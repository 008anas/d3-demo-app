from rest_framework import serializers

from .models import GeneticElement, Category


class GeneticElementSerializer(serializers.ModelSerializer):
    type = serializers.SerializerMethodField()

    class Meta:
        model = GeneticElement
        fields = ('id', 'type', 'glyph_thumbnail', 'default')

    def get_type(self, element):
        return element.name


class CategorySerializer(serializers.ModelSerializer):
    elements = GeneticElementSerializer(many=True)

    class Meta:
        model = Category
        fields = ('name', 'elements')