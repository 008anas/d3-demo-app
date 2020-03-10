from django.db import models

from app.genetic_element.models import GeneticElement
from app.specie.models import Specie


def species_dir(self, filename):
    # file will be uploaded to MEDIA_ROOT/matrices/<specie_name>/<filename>
    return 'matrices/' + self.specie.slug + '/' + filename


class Parameter(models.Model):
    name = models.CharField(max_length=255, help_text='Used for interface display')
    alias = models.CharField(max_length=255, unique=True, help_text='For internal usage. Must match checker list')
    specie = models.ForeignKey(Specie, on_delete=models.CASCADE)
    genetic_element = models.ForeignKey(GeneticElement, on_delete=models.CASCADE)
    matrix_file = models.FileField(null=True, blank=True, upload_to=species_dir, help_text='Matrix file')
    genome_min = models.FloatField(blank=True, default=0,
                                   help_text='If empty it will be automatically calculated from genome file')
    genome_max = models.FloatField(blank=True, default=1,
                                   help_text='If empty it will be automatically calculated from genome file')
    active = models.BooleanField(default=True, help_text='Determine if it can be used')
    created_at = models.DateTimeField('creation date', auto_now_add=True, editable=False)
    updated_at = models.DateTimeField('last update', auto_now=True, editable=False)

    class Meta:
        unique_together = ['alias', 'specie', 'genetic_element']

    def __str__(self):
        return self.name
