import sys

from Bio import SeqIO
from django.contrib import admin
from django.contrib.admin import ModelAdmin
from django.core.exceptions import ValidationError

from sqrutiny.settings import BASE_DIR
from .models import Parameter

sys.path.insert(0, BASE_DIR + '/../../dev')
from tools import checker


@admin.register(Parameter)
class GeneticElementAdmin(ModelAdmin):
    fieldsets = (
        ('Information',
         {'fields': (
             'name', 'alias', 'specie', 'genome_min', 'genome_max', 'genetic_element', 'description', 'matrix_file')}),
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

    def save_model(self, request, obj, form, change):
        if (not obj.genome_min or not obj.genome_max) and obj.specie.genome_gbk:
            try:
                handle = open(obj.specie.genome_gbk.path, 'rU')
                for record in SeqIO.parse(handle, 'genbank'):
                    result = checker(str(record.seq), parameter_dict={
                        obj.alias: dict(name=obj.name, min=0, max=1, elements={obj.alias},
                                        matrix='')})
                    if len(result) > 0:
                        raw = [x['raw_score'] for x in result[0]['scores']]
                        if not obj.genome_min:
                            obj.genome_min = min(raw)
                        if not obj.genome_max:
                            obj.genome_max = max(raw)
            except ValueError:
                raise ValidationError("Unable to calculate genome_min/ genome_max. Please enter it manually.")
        super().save_model(request, obj, form, change)
