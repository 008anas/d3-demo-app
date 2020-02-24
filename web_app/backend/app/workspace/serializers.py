import django_rq
from rest_framework import serializers

from .models import History
from ..construct.serializers import ConstructRetrieveSerializer


class HistorySerializer(serializers.ModelSerializer):
    id = serializers.SerializerMethodField()
    construct = ConstructRetrieveSerializer(read_only=True)
    status = serializers.SerializerMethodField()

    class Meta:
        model = History
        exclude = ['uuid', 'request_ip', 'deleted', 'updated_at']

    def get_url(self, history):
        request = self.context.get('request')
        return request.build_absolute_uri(str(history.uuid))

    def get_id(self, history):
        return history.uuid

    def get_status(self, obj):
        rq_job = django_rq.get_queue('default').fetch_job(str(obj.job_id))
        return rq_job.get_status()