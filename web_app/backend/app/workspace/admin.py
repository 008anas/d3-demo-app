from django.contrib import admin
from django.contrib.admin import ModelAdmin

from app.workspace.models import History


@admin.register(History)
class HistoryAdmin(ModelAdmin):
    fieldsets = (
        ('Information', {'fields': ('construct', 'job_id')}),
        ('Actions', {'fields': ('deleted',)})
    )
    add_fieldsets = (
        (None, {
            'classes': ('wide',),
            'fields': ('name',),
        }),
    )
    date_hierarchy = 'created_at'
    list_display = ('uuid', 'construct', 'job_id', 'deleted', 'created_at')
    search_fields = ('construct', 'job_id')
    ordering = ('created_at',)
