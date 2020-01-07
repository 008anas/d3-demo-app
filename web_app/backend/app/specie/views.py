from rest_framework import generics

from .models import Specie
from .serializers import SpecieSerializer


class SpecieListView(generics.ListAPIView):
    queryset = Specie.objects.filter(visible=True)
    serializer_class = SpecieSerializer


class SpecieRetrieveView(generics.RetrieveAPIView):
    lookup_field = 'slug'
    queryset = Specie.objects.filter(visible=True)
    serializer_class = SpecieSerializer
