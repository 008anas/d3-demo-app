from rest_framework import serializers

from .models import GeneticElement, Category


class GeneticElementSerializer(serializers.ModelSerializer):
    class Meta:
        model = GeneticElement
        fields = ('id', 'name', 'glyph_thumbnail', 'default')


class CategorySerializer(serializers.ModelSerializer):
    elements = GeneticElementSerializer(many=True)

    class Meta:
        model = Category
        fields = ('name', 'elements')