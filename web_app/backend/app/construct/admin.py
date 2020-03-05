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
        ('Information', {'fields': ('name', 'from_file', 'circular', 'dna_seq', 'protein_seq', 'specie', 'description')}),
        ('Actions', {'fields': ('deleted', 'example')})
    )
    add_fieldsets = (
        (None, {
            'classes': ('wide',),
            'fields': ('name',),
        }),
    )
    date_hierarchy = 'created_at'
    readonly_fields = ('from_file', 'protein_seq')
    list_display = ('uuid', 'name', 'dna_seq_length', 'protein_seq_length','tracks_count', 'circular', 'from_file', 'example', 'deleted', 'created_at')
    inlines = (TrackInline,)
    search_fields = ('name', 'sequence')
    ordering = ('-created_at',)

    def get_queryset(self, request):
        return Construct.objects.annotate(tracks_count=Count('tracks'))

    def tracks_count(self, obj):
        return obj.tracks_count

    def dna_seq_length(self, obj):
        return len(obj.dna_seq)

    def protein_seq_length(self, obj):
        return len(obj.protein_seq)
