from django.contrib import admin
from django.contrib.admin import ModelAdmin

from .models import Specie


@admin.register(Specie)
class SpecieAdmin(ModelAdmin):
    fieldsets = (
        ('Information', {'fields': ('name', 'tax_id', 'gc_content', 'slug', 'genome_gbk')}),
        ('Actions', {'fields': ('visible', 'default', 'comment')})
    )
    add_fieldsets = (
        (None, {
            'classes': ('wide',),
            'fields': ('name',),
        }),
    )
    list_display = ('name', 'slug', 'tax_id', 'default', 'created_at')
    search_fields = ('name', 'slug', 'tax_id')
    ordering = ('-created_at',)
