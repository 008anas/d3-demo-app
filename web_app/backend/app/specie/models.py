from django.db import models


class Specie(models.Model):
    name = models.CharField(max_length=255, unique=True)
    comment = models.TextField(blank=True, null=True, help_text='Short comment to display near name')
    slug = models.CharField(max_length=255, unique=True, help_text='Short name displayed in the url')
    ncbi_tax_id = models.IntegerField(unique=True)
    gc_content = models.FloatField(null=True, blank=True)
    codon_table = models.IntegerField(default=11)
    default = models.BooleanField(default=False, help_text='Determine which it will be used by default')
    visible = models.BooleanField(default=True, help_text='Determine if it is visible in the app')
    created_at = models.DateTimeField('creation date', auto_now_add=True, editable=False)
    updated_at = models.DateTimeField('last update', auto_now=True, editable=False)

    def __str__(self):
        return self.name

    def save(self, *args, **kwargs):
        if self.default:
            try:
                temp = Specie.objects.get(default=True)
                if self != temp:
                    temp.default = False
                    temp.save()
            except Specie.DoesNotExist:
                pass
        super(Specie, self).save(*args, **kwargs)
