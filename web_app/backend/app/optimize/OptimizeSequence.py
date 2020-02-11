import sys
import tempfile

import django_rq
# from sbol import *
from Bio import SeqIO
from rest_framework import status
from rest_framework.parsers import MultiPartParser
from rest_framework.response import Response
from rest_framework.views import APIView

from app.construct.models import Construct, Track
from app.construct.serializers import ConstructCreateSerializer
from app.genetic_element.models import GeneticElement
from app.matrix.models import Matrix
from app.serializers import GenBankSerializer
from app.specie.models import Specie
from app.utils import *
from app.workspace.models import History
from app.workspace.serializers import HistorySerializer
from sqrutiny.settings import BASE_DIR

sys.path.insert(0, BASE_DIR + '/../../dev')
from tools import checker, is_dna_seq_valid, match_sequence, seq_translator

TOOL_NAME = 'SQrutiny - Optimize Sequence - '
FEATURE_PREFIX = 'SQY_BOX_'


class OptimizeSequenceSkectherView(APIView):

    def post(self, request):
        serializer = ConstructCreateSerializer(data=request.data)

        if not serializer.is_valid(raise_exception=True):
            return Response(dict(msg='Invalid construct.'), status=status.HTTP_400_BAD_REQUEST)

        specie = Specie.objects.filter(tax_id=serializer.validated_data.get('specie_tax_id', None)).first()

        if specie is None:
            return Response({'msg': 'Specie with ncbi tax id ' + str(
                serializer.validated_data['specie_tax_id']) + ' was not found'}, status=status.HTTP_404_NOT_FOUND)

        construct = serializer.save(specie=specie)

        matrix = Matrix.objects.filter(active=True, specie=specie)

        if not matrix:
            return Response({'msg': 'Sorry but it was not possible to perform action. Please try later'},
                            status=status.HTTP_503_SERVICE_UNAVAILABLE)

        job = django_rq.enqueue(checker, sequence=construct.dna_seq, matrix_dict=self.matrix_to_dict(matrix),
                                circular=construct.circular, codon_table=specie.codon_table, result_ttl=-1)

        if job.get_status() == 'failed':
            return Response(dict(msg='Sorry There has been a problem. Please try later'),
                            status=status.HTTP_503_SERVICE_UNAVAILABLE)

        history = History.objects.create(
            name=TOOL_NAME + construct.name,
            construct=construct, job_id=job.id,
            request_ip=get_request_ip(request) or None
        )

        # Save history in current session
        self.request.session.setdefault('history', [])
        self.request.session['history'].append(str(history.uuid))
        self.request.session.modified = True

        serializer = HistorySerializer(instance=history, context={'request': request})

        return Response(serializer.data, status=status.HTTP_201_CREATED)

    @staticmethod
    def matrix_to_dict(matrix):
        return {entry.alias: entry.matrix_file.path for entry in matrix}


class OptimizeSequenceFileView(APIView):
    parser_classes = [MultiPartParser]

    def post(self, request):

        serializer = GenBankSerializer(data=request.data)

        serializer.is_valid(raise_exception=True)

        specie = Specie.objects.filter(tax_id=serializer.validated_data.get('specie_tax_id')).first()

        if specie is None:
            return Response({'msg': 'Specie with ncbi tax id ' + str(
                serializer.validated_data['specie_tax_id']) + ' was not found'},
                            status=status.HTTP_404_NOT_FOUND)

        file = serializer.validated_data['file']

        tmp_file = tempfile.NamedTemporaryFile()

        for chunk in file.chunks():
            tmp_file.write(chunk)
        tmp_file.flush()

        try:
            handle = open(tmp_file.name, 'rU')

            # try:
            for record in SeqIO.parse(handle, 'genbank'):

                if not record.seq:
                    return Response({'msg': 'No sequence was found. Please specify a valid sequence and try again.'},
                                    status=status.HTTP_400_BAD_REQUEST)

                tracks = []
                for feature in record.features:
                    if feature.type.lower() == 'cds':
                        tracks.append(dict(
                            genetic_element=GeneticElement.objects.filter(name__iexact=feature.type).first(),
                            sequence=str(feature.extract(record.seq)),
                            start=feature.location.nofuzzy_start,
                            end=feature.location.nofuzzy_end))

                construct = Construct.objects.create(
                    name=record.name,
                    dna_seq=record.seq,
                    protein_seq=seq_translator(str(record.seq)),
                    specie=specie,
                    description=record.annotations.get('description', None),
                    circular=True if record.annotations.get('topology', '').lower() == 'circular' else False,
                    from_file=True
                )

                if len(tracks):
                    Track.objects.bulk_create([Track(genetic_element=t.get('genetic_element'),
                                                     construct=construct,
                                                     sequence=t.get('sequence'),
                                                     start=t.get('start'),
                                                     end=t.get('end')) for t in tracks])

            # except:
            #     return Response({'msg': 'Invalid GenBank format file'},
            #                     status=status.HTTP_400_BAD_REQUEST)

        except:
            raise

        tmp_file.close()

        matrix = Matrix.objects.filter(active=True, specie=specie)

        if not matrix:
            return Response({'msg': 'Sorry but it was not possible to perform action. Please try later'},
                            status=status.HTTP_503_SERVICE_UNAVAILABLE)

        job = django_rq.enqueue(checker, sequence=construct.dna_seq, matrix_dict=self.matrix_to_dict(matrix),
                                circular=construct.circular, codon_table=specie.codon_table, result_ttl=-1)

        if job.get_status() == 'failed':
            return Response(dict(msg='Sorry There has been a problem. Please try later'),
                            status=status.HTTP_503_SERVICE_UNAVAILABLE)

        history = History.objects.create(
            name=TOOL_NAME + construct.name,
            construct=construct, job_id=job.id,
            request_ip=get_request_ip(request) or None
        )

        # Save history in current session
        self.request.session.setdefault('history', [])
        self.request.session['history'].append(str(history.uuid))
        self.request.session.modified = True

        serializer = HistorySerializer(instance=history, context={'request': request})

        return Response(serializer.data, status=status.HTTP_201_CREATED)

    @staticmethod
    def matrix_to_dict(matrix):
        return {entry.alias: entry.matrix_file.path for entry in matrix}


class SearchMotifView(APIView):

    def get(self, request):
        motif = self.request.query_params.get('motif')
        seq = self.request.query_params.get('sequence')
        if motif and seq:
            motif = motif.upper()
            if not is_dna_seq_valid(motif, False):
                return Response(dict(msg='Motif not valid. Please enter a valid one.'),
                                status=status.HTTP_400_BAD_REQUEST)
            seq = seq.upper()
            if not is_dna_seq_valid(seq, False):
                return Response(dict(msg='Sequence not valid. Please enter a valid one.'),
                                status=status.HTTP_400_BAD_REQUEST)

            match = match_sequence(motif, seq)
            return Response(dict(
                data=match,
                count=len(match)
            ), status=status.HTTP_200_OK)
        return Response(dict(msg='Invalid params'), status=status.HTTP_400_BAD_REQUEST)
