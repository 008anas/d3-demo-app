from rest_framework import generics

from .models import Category, GeneticElement
from .serializers import GeneticElementSerializer, CategorySerializer


class GeneticElementListView(generics.ListAPIView):
    queryset = GeneticElement.objects.filter(visible=True)
    serializer_class = GeneticElementSerializer


class CategoryListView(generics.ListAPIView):
    queryset = Category.objects.all()
    serializer_class = CategorySerializer
