from django.contrib import admin
from django.contrib.admin import ModelAdmin

from app.contact.models import Contact


@admin.register(Contact)
class SpecieAdmin(ModelAdmin):
    fieldsets = (
        ('Information', {'fields': ('name', 'email', 'subject', 'message')}),
    )
    add_fieldsets = (
        (None, {
            'classes': ('wide',),
            'fields': ('name',),
        }),
    )
    list_display = ('name', 'email', 'subject', 'message', 'created_at')
    search_fields = ('name', 'email', 'subject', 'message')
    ordering = ('created_at',)