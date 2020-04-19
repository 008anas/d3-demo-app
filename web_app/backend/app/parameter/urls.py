from django.urls import path

from app.parameter.views import ParametersListView

urlpatterns = [
    path('', ParametersListView.as_view())
]