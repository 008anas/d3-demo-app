from django.contrib import admin
from django.contrib.admin import ModelAdmin

from app.matrix.models import Matrix


@admin.register(Matrix)
class GeneticElementAdmin(ModelAdmin):
    fieldsets = (
        ('Information', {'fields': ('name', 'alias', 'specie', 'genetic_element', 'matrix_file')}),
        ('Actions', {'fields': ('active',)})
    )
    add_fieldsets = (
        (None, {
            'classes': ('wide',),
            'fields': ('name',),
        }),
    )
    list_display = ('name', 'specie', 'genetic_element', 'active', 'created_at')
    search_fields = ('name',)
    ordering = ('-created_at',)