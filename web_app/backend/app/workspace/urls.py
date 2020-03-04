from django.urls import path, include
from rest_framework import routers

from app.workspace.views import HistoryViewSet, HistoryCountView, PollJobView, ClearHistoryView, HistoryRetrieveView, \
    ExportThresholdView

router = routers.DefaultRouter(trailing_slash=False)
router.register(r'', HistoryViewSet, basename='History')


urlpatterns = [
    path('all', ClearHistoryView.as_view(), name='Clear user history'),
    path('count', HistoryCountView.as_view(), name="User workspace history count"),
    path('job/<uuid:job_id>', PollJobView.as_view()),
    path('history/<uuid:history_id>', HistoryRetrieveView.as_view()), # Retrieve history without throwing 404 error
    path('export/<uuid:history_id>', ExportThresholdView.as_view()),  # Retrieve history without throwing 404 error
    path('', include(router.urls))
]
