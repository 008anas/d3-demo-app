from rest_framework import serializers

from app.construct.serializers import ConstructRetrieveSerializer
from .models import History


class HistorySerializer(serializers.ModelSerializer):
    id = serializers.SerializerMethodField()
    construct = ConstructRetrieveSerializer(read_only=True)

    class Meta:
        model = History
        exclude = ['uuid', 'request_ip', 'deleted', 'updated_at']

    def get_url(self, history):
        request = self.context.get('request')
        return request.build_absolute_uri(str(history.uuid))

    def get_id(self, history):
        return history.uuid