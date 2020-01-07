from django.contrib import admin
from django.contrib.admin import ModelAdmin
from django.db.models import Count

from .models import Construct, Track


class TrackInline(admin.TabularInline):
    model = Track
    extra = 1


@admin.register(Construct)
class ConstructAdmin(ModelAdmin):
    fieldsets = (
        ('Information', {'fields': ('label', 'sequence', 'specie')}),
        ('Actions', {'fields': ('deleted', 'example')})
    )
    add_fieldsets = (
        (None, {
            'classes': ('wide',),
            'fields': ('name',),
        }),
    )
    date_hierarchy = 'created_at'
    list_display = ('label', 'sequence_length', 'tracks_count', 'deleted', 'created_at')
    inlines = (TrackInline,)
    search_fields = ('label', 'sequence')
    ordering = ('created_at',)

    def get_queryset(self, request):
        return Construct.objects.annotate(tracks_count=Count('tracks'))

    def tracks_count(self, obj):
        return obj.tracks_count

    def sequence_length(self, obj):
        return len(obj.sequence)
