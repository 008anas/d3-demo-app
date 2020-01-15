from rest_framework import viewsets, serializers, status, generics
from rest_framework.response import Response
from rest_framework.views import APIView

from app.construct.models import Construct
from app.construct.serializers import ConstructRetrieveSerializer, ConstructCreateSerializer


class ConstructListRetrieveView(viewsets.ModelViewSet):
    queryset = Construct.objects.filter(deleted=False).order_by('-created_at')
    serializer_class = ConstructRetrieveSerializer
    lookup_field = 'uuid'


class ConstructExampleView(generics.ListCreateAPIView):
    queryset = Construct.objects.filter(deleted=False, example=True)
    serializer_class = ConstructRetrieveSerializer


class SerializerClass(serializers.Serializer):
    option = serializers.CharField()
    construct = ConstructCreateSerializer()


class ExportConstructView(APIView):

    def post(self, request):
        serializer = SerializerClass(data=request.data)

        if not serializer.is_valid():
            return Response(dict(msg='Invalid construct.'), status=status.HTTP_400_BAD_REQUEST)

        construct = serializer.validated_data.get('construct')



        return Response(serializer.data, status=status.HTTP_201_CREATED)
