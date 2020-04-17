from django.urls import path, include
from rest_framework import routers

from app.parameter.views import ParametersListView

urlpatterns = [
    path('', ParametersListView.as_view())
]