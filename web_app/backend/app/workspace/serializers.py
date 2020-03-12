import django_rq
from rest_framework import serializers

from app.parameter.models import Parameter
from .models import History
from ..construct.serializers import ConstructRetrieveSerializer


class HistorySerializer(serializers.ModelSerializer):
    id = serializers.SerializerMethodField()
    construct = ConstructRetrieveSerializer(read_only=True)
    status = serializers.SerializerMethodField()

    class Meta:
        model = History
        exclude = ['uuid', 'deleted', 'updated_at']

    def get_url(self, history):
        request = self.context.get('request')
        return request.build_absolute_uri(str(history.uuid))

    def get_id(self, history):
        return history.uuid

    def get_status(self, obj):
        rq_job = django_rq.get_queue('default').fetch_job(str(obj.job_id))
        return rq_job.get_status()


class FilterSerializer(serializers.Serializer):
    key = serializers.ChoiceField(choices=set(Parameter.objects.filter(active=True).values_list('alias', flat=True)))
    value = serializers.FloatField()
    op = serializers.ChoiceField(choices=['<', '<=', '=', '>', '>='])
    type = serializers.ChoiceField(choices=['raw', 'norm'])


class ExportResultsSerializer(serializers.Serializer):
    filters = FilterSerializer(many=True, required=False)
    bulk = serializers.BooleanField(required=False)
