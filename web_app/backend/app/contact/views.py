from django.conf import settings
from django.core.mail import send_mail
from rest_framework import generics, status
from rest_framework.response import Response

from .serializers import ContactSerializer


class ContactCreateView(generics.CreateAPIView):

    def post(self, request):
        comment = ContactSerializer(data=request.data)

        if not comment.is_valid():
            return Response({'serializer': comment})

        form_subject = comment.data['subject']
        form_email = comment.data['email']
        form_message = comment.data['message']
        form_name = comment.data['name']

        send_mail(form_subject,
                  form_message,
                  form_email,
                  [settings.EMAIL_HOST_USER],
                  fail_silently=False
                  )
        return Response(comment.data, status=status.HTTP_201_CREATED)