"""SQRUTINY URL Configuration"""

from django.contrib import admin
from django.urls import path, include
from django.conf.urls.static import static

from . import settings


urlpatterns = [
    path('admin/', admin.site.urls),
    path('admin/django-rq', include('django_rq.urls')),
    path('api/v1/', include('app.urls'))
]

if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL,
                          document_root=settings.MEDIA_ROOT)
