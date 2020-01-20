from rest_framework import viewsets, generics

from app.construct.models import Construct
from app.construct.serializers import ConstructRetrieveSerializer


class ConstructListRetrieveView(viewsets.ModelViewSet):
    queryset = Construct.objects.filter(deleted=False).order_by('-created_at')
    serializer_class = ConstructRetrieveSerializer
    lookup_field = 'uuid'


class ConstructExampleView(generics.ListCreateAPIView):
    queryset = Construct.objects.filter(deleted=False, example=True)
    serializer_class = ConstructRetrieveSerializer
