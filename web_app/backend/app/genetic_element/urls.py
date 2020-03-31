from django.urls import path

from .views import GeneticElementListView, CategoryListView

urlpatterns = [
    path('', GeneticElementListView.as_view(), name='Tracks'),
    path('categories', CategoryListView.as_view(), name='Categories with tracks')
]
