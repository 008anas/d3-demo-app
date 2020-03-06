import operator

import django_rq
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from rest_framework import viewsets, status, serializers
from rest_framework.response import Response
from rest_framework.views import APIView

from .models import History
from .serializers import HistorySerializer


class HistoryCountView(APIView):

    def get(self, request):
        return Response(dict(count=len(self.request.session.get('history', []))))


class HistoryRetrieveView(APIView):

    def get(self, request, history_id):
        if str(history_id) not in request.session.get('history', []):
            return Response(dict(msg='History with such id not found.'), status=status.HTTP_422_UNPROCESSABLE_ENTITY)

        history = History.objects.filter(uuid=str(history_id), deleted=False).first()

        serializer = HistorySerializer(instance=history, context={'request': request})

        return Response(serializer.data, status=status.HTTP_200_OK)


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
        return History.objects.filter(deleted=False).order_by('-created_at')
        # return History.objects.filter(uuid__in=self.request.session.get('history', []), deleted=False).order_by(
        #   '-created_at')


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

        return Response(dict(
            status=rq_job.get_status(),
            result=rq_job.result
        ), status=status.HTTP_200_OK)


class ThresholdSerializer(serializers.Serializer):
    start = serializers.IntegerField()
    end = serializers.IntegerField()
    raw_score = serializers.FloatField()
    norm_score = serializers.FloatField()


class ExportThresholdView(APIView):

    def post(self, request, history_id):

        op_fn = {
            '<': operator.lt,
            '<=': operator.le,
            '=': operator.eq,
            '>': operator.gt,
            '>=': operator.ge,
        }

        scores_type = ['raw', 'norm']
        modes = ['default', 'bulk']

        # if str(history_id) not in request.session.get('history', []):
        #   return Response(dict(msg='History with such id not found.'), status=status.HTTP_422_UNPROCESSABLE_ENTITY)

        history = History.objects.filter(uuid=str(history_id), deleted=False).first()

        job = django_rq.get_queue('default').fetch_job(str(history.job_id))

        record = history.construct.to_record()

        if not record or not job:
            return Response(dict(msg='Unable to export'), status=status.HTTP_422_UNPROCESSABLE_ENTITY)

        result = job.result

        filters = request.data.get('filters', [])
        mode = request.data.get('mode', 'default')

        for f in filters:
            if f.get('key', None) in result and f.get('value') is not None and f.get('op', None) in op_fn and \
                    f.get('type') in scores_type:
                # Override results
                result[f['key']] = [score for score in result[f['key']] if
                                    op_fn.get(f['op'])(
                                        score.get('raw_score') if f['type'] == 'raw' else score.get('norm_score'),
                                        float(f['value']))]
        # Remove non filters
        if mode == 'default':
            result = {k: result[k] for k in set([d.get('key', None) for d in filters])}

        # Append result scores
        for r in result:
            for t in result[r]:
                feature = SeqFeature(FeatureLocation(t.get('start') - 1, t.get('end')), type='SQY_SCORE')
                feature.qualifiers = dict(norm_score=t.get('norm_score'), raw_score=t.get('raw_score'))
                feature.qualifiers['sqy_type'] = r
                record.features.append(feature)

        # Save as GenBank file
        file = open('example.gb', 'w')
        SeqIO.write(record, file, 'genbank')

        return Response(dict(msg='GenBank generated'), status=status.HTTP_200_OK)
