import uuid

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from django.db import models

from app.genetic_element.models import GeneticElement
from app.specie.models import Specie


class Construct(models.Model):
    uuid = models.UUIDField(editable=False, unique=True, default=uuid.uuid4)
    name = models.CharField(max_length=255)
    dna_seq = models.TextField(blank=True)
    protein_seq = models.TextField(blank=True)
    tracks = models.ManyToManyField(GeneticElement, through='Track')
    specie = models.ForeignKey(Specie, on_delete=models.CASCADE)
    circular = models.BooleanField(blank=True, default=False)
    description = models.TextField(blank=True, null=True)
    example = models.BooleanField(default=False, help_text='Determine if is it used as example')
    from_file = models.BooleanField(default=False,
                                    help_text='Describe if construct was loaded from file or using SQRUTINY sketcher')
    deleted = models.BooleanField(default=False, help_text='Determine if it is visible in the application')
    created_at = models.DateTimeField('creation date', auto_now_add=True)
    updated_at = models.DateTimeField('last update', auto_now=True)

    def __str__(self):
        return self.name or ''

    def save(self, *args, **kwargs):
        if self.example:
            try:
                temp = Construct.objects.get(example=True)
                if self != temp:
                    temp.default = False
                    temp.save()
            except Construct.DoesNotExist:
                pass
        super(Construct, self).save(*args, **kwargs)

    def from_record(self, record):
        return

    def to_record(self):
        # Create sequence
        sequence_object = Seq(self.dna_seq, IUPAC.unambiguous_dna)

        # Create record
        record = SeqRecord(sequence_object,
                           name='Example_GenBank',
                           description=self.description or '')

        record.annotations['organism'] = record.annotations['source'] = self.specie.name
        record.annotations['topology'] = 'circular' if self.circular else 'linear'

        first_ft = SeqFeature(FeatureLocation(0, len(self.dna_seq)), type='source')
        first_ft.qualifiers['organism'] = self.specie.name
        first_ft.qualifiers['db_xref'] = 'taxon:' + str(self.specie.tax_id)

        record.features.append(first_ft)

        for t in self.track_set.all():
            feature = SeqFeature(FeatureLocation(t.start-1, t.end), type='SQY_BOX')
            feature.qualifiers['sqy_type'] = t.genetic_element.name
            if t.label:
                feature.qualifiers['label'] = t.label
            record.features.append(feature)

        return record


class Track(models.Model):
    genetic_element = models.ForeignKey(GeneticElement, on_delete=models.CASCADE)
    construct = models.ForeignKey(Construct, on_delete=models.CASCADE)
    label = models.CharField(max_length=255, null=True, blank=True)
    sequence = models.TextField(help_text='DNA sequence')
    start = models.IntegerField(blank=True)
    end = models.IntegerField(blank=True)
    color = models.CharField(max_length=10, default='#4e0a77')

    def __str__(self):
        return self.label or ''
