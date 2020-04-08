import operator
import tempfile

import django_rq
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from rest_framework import viewsets, status
from rest_framework.response import Response
from rest_framework.views import APIView

from .models import History
from .serializers import HistorySerializer, ExportResultsSerializer


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
        return History.objects.filter(uuid__in=self.request.session.get('history', []), deleted=False).order_by(
            '-created_at')


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


class ExportThresholdView(APIView):

    def post(self, request, history_id):

        serializer = ExportResultsSerializer(data=request.data)

        serializer.is_valid(raise_exception=True)

        op_fn = {
            '<': operator.lt,
            '<=': operator.le,
            '=': operator.eq,
            '>': operator.gt,
            '>=': operator.ge
        }

        if str(history_id) not in request.session.get('history', []):
            return Response(dict(msg='History with such id not found.'), status=status.HTTP_422_UNPROCESSABLE_ENTITY)

        history = History.objects.filter(uuid=str(history_id), deleted=False).first()

        job = django_rq.get_queue('default').fetch_job(str(history.job_id))

        record = history.construct.to_record()

        if not record or not job:
            return Response(dict(msg='Unable to export'), status=status.HTTP_422_UNPROCESSABLE_ENTITY)

        result = job.result

        filters = serializer.validated_data.get('filters', [])
        bulk = serializer.validated_data.get('bulk', False)

        keys = [f.get('key', '') for f in filters]

        # Keep only in filters (if bulk mode)
        if not bulk:
            result = [r for r in result if r['alias'] in keys]

        for f in filters:
            # Override results
            single = next((r for r in result if r['alias'] == f['key']), None)
            if single:
                single['scores'] = [s for s in single['scores'] if op_fn.get(f['op'])(
                    s.get('raw_score') if f['type'] == 'raw' else s.get('norm_score'), float(f['value']))]
                result[result.index(single)] = single

        # Append result scores
        for r in result:
            for t in r.get('scores', []):
                feature = SeqFeature(FeatureLocation(t['start'], t['end']), type='SQY_SCORE')
                feature.qualifiers = dict(norm_score=t['norm_score'], raw_score=t['raw_score'])
                feature.qualifiers['sqy_type'] = r['alias']
                record.features.append(feature)

        tmp_file = tempfile.TemporaryFile(mode="r+")

        SeqIO.write(record, tmp_file, 'genbank')

        tmp_file.seek(0)
        data = tmp_file.read()
        tmp_file.close()

        file = {}

        if data:
            file = dict(
                data=data,
                filename=history.name + '_scores_export.gb',
                mimetype='application/genbank')

        return Response(file, status=status.HTTP_200_OK)
