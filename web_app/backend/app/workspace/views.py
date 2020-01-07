import django_rq
from rest_framework import viewsets, status
from rest_framework.response import Response
from rest_framework.views import APIView

from .models import History
from .serializers import HistorySerializer


class HistoryCountView(APIView):

    def get(self, request):
        return Response(dict(count=len(self.request.session.get('history'))))


class HistoryViewSet(viewsets.ModelViewSet):
    serializer_class = HistorySerializer
    lookup_field = 'uuid'

    def destroy(self, request, *args, **kwargs):

        history = self.get_object()
        history.deleted = True
        history.save()
        self.request.session['history'].remove(kwargs.get('uuid'))
        self.request.session.modified = True

        return Response(dict(msg='History deleted successfully!'), status=status.HTTP_204_NO_CONTENT)

    def get_queryset(self):
        return History.objects.filter(uuid__in=self.request.session.get('history', []), deleted=False)


class ClearHistoryView(APIView):

    def delete(self, request):

        self.get_queryset().update(deleted=True)

        self.request.session.clear()

        return Response(dict(msg='History cleared successfully!'), status=status.HTTP_200_OK)

    def get_queryset(self):
        return History.objects.filter(uuid__in=self.request.session.get('history', []), deleted=False)


class PollJobView(APIView):

    def get(self, request, job_id):
        rq_job = django_rq.get_queue('default').fetch_job(str(job_id))

        if rq_job is None:
            return Response(dict(msg='Unable to find job'), status=status.HTTP_404_NOT_FOUND)

        job = dict(status=rq_job.get_status(), result=rq_job.result, progress=rq_job.meta.get('progress', None))

        return Response(job, status=status.HTTP_200_OK)
