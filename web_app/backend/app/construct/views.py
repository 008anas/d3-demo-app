import tempfile

from Bio import SeqIO
from rest_framework import viewsets, generics, status
from rest_framework.parsers import MultiPartParser
from rest_framework.response import Response
from rest_framework.views import APIView

from app.construct.models import Construct
from app.construct.serializers import ConstructRetrieveSerializer
from app.serializers import GenBankSerializer

FEATURE_PREFIX = 'SQY_BOX'
CDS_COLOR = '#4e0a77'
FEATURES_COLOR = '#f0ca68'
DUMMY_COLOR = '#62382f'


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

            try:
                for record in SeqIO.parse(handle, 'genbank'):

                    if not record.seq or not len(record.seq):
                        return Response(
                            dict(msg='No sequence was found. Please specify a valid sequence and try again.'),
                            status=status.HTTP_400_BAD_REQUEST
                        )

                    if not record.features or not len(record.features):
                        return Response(
                            dict(msg='No feature was found. Please specify some CDS features and try again.'),
                            status=status.HTTP_400_BAD_REQUEST
                        )

                    tracks = []
                    sqy_tracks = []

                    if record.features[0].type.lower() == 'source':
                        record.features.pop(0)

                    tracks.append(dict(
                        label=record.features[0].qualifiers.get('locus_tag', ['Track 1'])[0],
                        type='CDS' if record.features[0].type.lower() == 'cds' else 'Dummy',
                        sequence=str(record.features[0].extract(record.seq)),
                        color=CDS_COLOR,
                        start=record.features[0].location.nofuzzy_start,
                        end=record.features[0].location.nofuzzy_end
                    ))

                    record.features.pop(0)

                    i = 1
                    for feature in record.features:
                        last_item = tracks[-1]
                        if feature.type.lower() == 'cds':
                            if last_item['end'] > feature.location.nofuzzy_start:  # No overlapping
                                return Response(dict(
                                    msg='The file contains overlapped features. '
                                        'Overlapping is not supported. '
                                        'Please correct the file and upload again.'),
                                    status=status.HTTP_400_BAD_REQUEST)
                            tracks.append(dict(
                                label=feature.qualifiers.get('locus_tag', 'Track ' + str(len(tracks)+1)),
                                type=feature.type,
                                sequence=str(feature.extract(record.seq)),
                                color=CDS_COLOR,
                                start=feature.location.nofuzzy_start,
                                end=feature.location.nofuzzy_end))
                        elif feature.type.upper() == FEATURE_PREFIX:
                            sqy_tracks.append(dict(
                                label=feature.qualifiers.get('locus_tag', 'Track ' + str(len(tracks)+1)),
                                type=feature.type,
                                sequence=str(feature.extract(record.seq)),
                                color=FEATURES_COLOR,
                                start=feature.location.nofuzzy_start,
                                end=feature.location.nofuzzy_end))
                        else:
                            if last_item.get('type', '').lower() != 'dummy':
                                tracks.append(dict(
                                    label=feature.qualifiers.get('locus_tag', 'Track ' + str(len(tracks)+1)),
                                    type='Dummy',
                                    sequence=str(feature.extract(record.seq)),
                                    color=DUMMY_COLOR,
                                    start=last_item['end'] + 1,
                                    end=feature.location.nofuzzy_end))
                            else:
                                tracks[-1]['end'] = feature.location.nofuzzy_end
                        i += 1

                    data = dict(
                        name=record.name,
                        dna_seq=str(record.seq),
                        description=record.annotations.get('description', None),
                        circular=True if record.annotations.get('topology', '').lower() == 'circular' else False,
                        from_file=True
                    )

                    if len(tracks):
                        data['tracks'] = tracks

                    if len(sqy_tracks):
                        data['sqy_tracks'] = sqy_tracks

            except:
                return Response(dict(msg='Invalid GenBank format file'),
                                status=status.HTTP_400_BAD_REQUEST)

        except:
            return Response(dict(msg='Error while handling. Please try with different file.'),
                            status=status.HTTP_400_BAD_REQUEST)

        tmp_file.close()

        return Response(data, status=status.HTTP_200_OK)
