from django.urls import path

from .views import SpecieListView, SpecieRetrieveView

urlpatterns = [
    path('', SpecieListView.as_view()),
    path('/<slug>', SpecieRetrieveView.as_view())
]
