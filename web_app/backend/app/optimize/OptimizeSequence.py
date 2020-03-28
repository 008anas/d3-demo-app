import sys

import django_rq
from rest_framework import status
from rest_framework.response import Response
from rest_framework.views import APIView
from rest_framework_tracking.mixins import LoggingMixin

from app.construct.serializers import ConstructCreateSerializer
from app.parameter.models import Parameter
from app.specie.models import Specie
from app.workspace.models import History
from app.workspace.serializers import HistorySerializer
from sqrutiny.settings import BASE_DIR

sys.path.insert(0, BASE_DIR + '/../../dev')
from tools import checker, is_dna_seq_valid, match_sequence


TOOL_NAME = 'SQrutiny - '


class OptimizeSequenceSkectherView(LoggingMixin, APIView):

    def post(self, request):
        serializer = ConstructCreateSerializer(data=request.data)

        if not serializer.is_valid(raise_exception=True):
            return Response(dict(msg='Invalid construct.'), status=status.HTTP_400_BAD_REQUEST)

        tax_id = serializer.validated_data.get('specie_tax_id', None)

        specie = Specie.objects.filter(tax_id=tax_id).first()

        if specie is None:
            return Response(dict(msg='Specie with ncbi tax id ' + str(tax_id) + ' was not found'),
                            status=status.HTTP_404_NOT_FOUND)

        construct = serializer.save(specie=specie)

        parameters = Parameter.objects.filter(active=True, specie=specie).values_list('id', flat=True)

        if not parameters:
            return Response(dict(msg='Sorry but it was not possible to perform action. Please try later'),
                            status=status.HTTP_503_SERVICE_UNAVAILABLE)

        job = django_rq.enqueue(checker, sequence=construct.dna_seq,
                                parameter_dict=Parameter.to_dict_by_id(ids=parameters),
                                circular=construct.circular,
                                codon_table=specie.codon_table, result_ttl=-1)

        if job.get_status() == 'failed':
            return Response(dict(msg='Sorry There has been a problem. Please try later'),
                            status=status.HTTP_503_SERVICE_UNAVAILABLE)

        history = History.objects.create(
            name=TOOL_NAME + construct.name,
            construct=construct, job_id=job.id
        )

        # save the current history in session
        request.session.setdefault('history', [])
        request.session['history'].append(str(history.uuid))
        request.session.modified = True

        serializer = HistorySerializer(instance=history, context={'request': request})

        return Response(serializer.data, status=status.HTTP_201_CREATED)


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
