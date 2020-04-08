from django.conf import settings
from django.core.mail import send_mail
from rest_framework import generics, status
from rest_framework.response import Response

from .serializers import ContactSerializer


class ContactCreateView(generics.CreateAPIView):

    def post(self, request):
        contact = ContactSerializer(data=request.data)

        contact.is_valid(raise_exception=True)

        send_mail(contact.data.get('subject', ''),
                  contact.data.get('message'),
                  contact.data.get('email'),
                  [settings.EMAIL_HOST_USER],
                  fail_silently=False
                  )
        return Response(dict(msg='Ok'), status=status.HTTP_201_CREATED)
