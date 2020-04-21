from rest_framework import serializers

from app.parameter.models import Parameter


class ParameterSerializer(serializers.ModelSerializer):
    class Meta:
        model = Parameter
        fields = ['name', 'alias', 'description', 'genome_min', 'genome_max']