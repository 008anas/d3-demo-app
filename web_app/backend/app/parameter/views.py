from rest_framework import generics

from app.parameter.models import Parameter
from app.parameter.serializers import ParameterSerializer


class ParametersListView(generics.ListCreateAPIView):
    serializer_class = ParameterSerializer

    def get_queryset(self):
        queryset = Parameter.objects.filter(active=True)
        specie_id = self.request.query_params.get('specie_id', None)
        if specie_id is not None:
            queryset = queryset.filter(specie__tax_id=specie_id)
        return queryset
