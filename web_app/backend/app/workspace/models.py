import uuid

from django.db import models

from app.construct.models import Construct


class History(models.Model):
    uuid = models.UUIDField(editable=False, unique=True, default=uuid.uuid4)
    name = models.CharField(max_length=255, blank=True)
    construct = models.ForeignKey(Construct, on_delete=models.CASCADE)
    job_id = models.UUIDField()
    request_ip = models.GenericIPAddressField(null=True)
    deleted = models.BooleanField(default=False)
    created_at = models.DateTimeField('creation date', auto_now_add=True)
    updated_at = models.DateTimeField('last update', auto_now=True)

    class Meta:
        verbose_name_plural = 'Histories'