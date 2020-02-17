from django.core.validators import FileExtensionValidator
from rest_framework import serializers


class GenBankSerializer(serializers.Serializer):
    file = serializers.FileField(validators=[FileExtensionValidator(allowed_extensions=['genbank', 'gbk', 'gb'])])