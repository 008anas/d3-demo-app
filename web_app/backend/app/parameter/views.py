from rest_framework import generics

from app.parameter.models import Parameter
from app.parameter.serializers import ParameterSerializer


class ParametersListView(generics.ListCreateAPIView):
    queryset = Parameter.objects.filter(active=True)
    serializer_class = ParameterSerializer