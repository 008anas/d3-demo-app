import tempfile

from Bio import SeqIO
from rest_framework import viewsets, generics, status
from rest_framework.parsers import MultiPartParser
from rest_framework.response import Response
from rest_framework.views import APIView

from app.construct.models import Construct
from app.construct.serializers import ConstructRetrieveSerializer
from app.serializers import GenBankSerializer

FEATURE_PREFIX = 'SQY_BOX_'


class ConstructListRetrieveView(viewsets.ModelViewSet):
    queryset = Construct.objects.filter(deleted=False).order_by('-created_at')
    serializer_class = ConstructRetrieveSerializer
    lookup_field = 'uuid'


class ConstructExampleView(generics.ListCreateAPIView):
    queryset = Construct.objects.filter(deleted=False, example=True)
    serializer_class = ConstructRetrieveSerializer


class FromGenBankView(APIView):
    parser_classes = [MultiPartParser]

    def post(self, request):

        serializer = GenBankSerializer(data=request.data)

        serializer.is_valid(raise_exception=True)

        file = serializer.validated_data['file']

        tmp_file = tempfile.NamedTemporaryFile()

        for chunk in file.chunks():
            tmp_file.write(chunk)
        tmp_file.flush()

        try:
            handle = open(tmp_file.name, 'rU')

            # try:
            for record in SeqIO.parse(handle, 'genbank'):

                if not record.seq or not len(record.seq):
                    return Response({'msg': 'No sequence was found. Please specify a valid sequence and try again.'},
                                    status=status.HTTP_400_BAD_REQUEST)

                tracks = []
                sqy_tracks = []

                for feature in record.features:
                    if feature.type.lower() == 'cds':
                        tracks.append(dict(
                            type=feature.type,
                            sequence=str(feature.extract(record.seq)),
                            start=feature.location.nofuzzy_start,
                            end=feature.location.nofuzzy_end))
                    elif feature.type.upper().startswith(FEATURE_PREFIX):
                        sqy_tracks.append(dict(
                            type=feature.type,
                            sequence=str(feature.extract(record.seq)),
                            start=feature.location.nofuzzy_start,
                            end=feature.location.nofuzzy_end))

                data = dict(
                    name=record.name,
                    dna_seq=str(record.seq),
                    description=record.annotations.get('description', None),
                    circular=True if record.annotations.get('topology', '').lower() == 'circular' else False
                )

                if len(tracks):
                    data['tracks'] = tracks

            # except:
            #     return Response({'msg': 'Invalid GenBank format file'},
            #                     status=status.HTTP_400_BAD_REQUEST)

        except:
            raise

        tmp_file.close()

        return Response(data, status=status.HTTP_201_CREATED)
