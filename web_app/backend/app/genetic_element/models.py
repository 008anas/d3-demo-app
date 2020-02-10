import os
import xml.etree.cElementTree as et

from django.core.exceptions import ValidationError
from django.db import models
from django.utils.safestring import mark_safe
from django.utils.text import slugify
from django.utils.translation import ugettext_lazy as _

from app.specie.models import Specie


class Category(models.Model):
    name = models.CharField(max_length=255, unique=True)
    created_at = models.DateTimeField('creation date', auto_now_add=True)
    updated_at = models.DateTimeField('last update', auto_now=True)

    class Meta:
        verbose_name = 'Category'
        verbose_name_plural = 'Categories'
        ordering = ['id']

    def __str__(self):
        return self.name


@mark_safe
def image_preview(obj):
    if obj.pk:  # if object has already been saved and has a primary key, show picture preview
        return '<a href="{src}" target="_blank"><img src="{src}" alt="{title}" width={width} height={height} /></a>'.format(
            width=250,
            height=200,
            src=obj.image.url,
            title=obj.name,
        )
    return _('(choose a valid image, save and continue editing to see the preview)')


image_preview.allow_tags = True
image_preview.short_description = _('Image Preview')


@mark_safe
def image_list_preview(obj):
    return '<a href="{src}" target="_blank"><img src="{src}" alt="{title}" width={width} height={height} /></a>'.format(
        width=150,
        height=100,
        src=obj.glyph_thumbnail.url,
        title=obj.name
    )


image_preview.allow_tags = True
image_preview.short_description = _('Image Preview')


def upload_image(instance, filename):
    filename_base, filename_ext = os.path.splitext(filename)
    return 'elem_img/{filename}{extension}'.format(
        filename=slugify(filename_base),
        extension=filename_ext.lower(),
    )


def validate_svg(f):
    # Find 'start' word in file and get 'tag' from there
    f.seek(0)
    tag = None
    try:
        for event, el in et.iterparse(f, ('start',)):
            tag = el.tag
            break
    except et.ParseError:
        pass

    # Check that this 'tag' is correct
    if tag != '{http://www.w3.org/2000/svg}svg':
        raise ValidationError('Uploaded file is not a valid SVG file.')

    # 'reset' file
    f.seek(0)

    return f


class GeneticElement(models.Model):
    name = models.CharField(max_length=255, unique=True)
    glyph_thumbnail = models.FileField('Genetic Element Icon', upload_to=upload_image, validators=[validate_svg])
    category = models.ForeignKey(Category, models.SET_NULL, related_name='elements', blank=True, null=True)
    role = models.URLField(unique=True, null=True, help_text='Uniform Resource Identifier (URI). Found here: http://www.sequenceontology.org/browser/obob.cgi')
    visible = models.BooleanField(default=True, help_text='Determine if it is visible in the application')
    default = models.BooleanField(default=False, help_text='Determine if it will be used by default')
    order = models.IntegerField(unique=True)
    created_at = models.DateTimeField('creation date', auto_now_add=True)
    updated_at = models.DateTimeField('last update', auto_now=True)

    class Meta:
        verbose_name = 'Genetic Element'
        verbose_name_plural = 'Genetic Elements'

    def __str__(self):
        return self.name or ''