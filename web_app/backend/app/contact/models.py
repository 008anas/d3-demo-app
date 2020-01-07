from django.db import models


class Contact(models.Model):
    name = models.CharField(max_length=255, blank=True, help_text='Sender name')
    email = models.EmailField(help_text='Sender email')
    subject = models.CharField(blank=True, max_length=255, help_text='Email subject')
    message = models.TextField(help_text='Email body')
    created_at = models.DateTimeField('creation date', auto_now_add=True)
    updated_at = models.DateTimeField('last update', auto_now=True)
