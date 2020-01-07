from django.urls import path, include
from rest_framework import routers

from .views import ConstructListRetrieveView, ExportConstructView, ConstructExampleView

router = routers.DefaultRouter()
router.register(r'', ConstructListRetrieveView)

urlpatterns = [
    path('export', ExportConstructView.as_view()),
    path('example', ConstructExampleView.as_view()),
    path('', include(router.urls))
]
