from rest_framework import serializers

from .models import Specie


class SpecieSerializer(serializers.ModelSerializer):
    class Meta:
        model = Specie
        fields = ('name', 'comment','slug', 'tax_id', 'tax_link', 'default')
