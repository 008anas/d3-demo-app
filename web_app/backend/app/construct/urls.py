from django.urls import path, include
from rest_framework import routers

from .views import ConstructListRetrieveView, ConstructExampleView

router = routers.DefaultRouter()
router.register(r'', ConstructListRetrieveView, basename='Construct')

urlpatterns = [
    path('example', ConstructExampleView.as_view()),
    path('', include(router.urls))
]
