# Generated by Django 2.2.8 on 2020-02-10 09:52

import django.db.models.deletion
from django.db import migrations, models

import app.genetic_element.models


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Category',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=255, unique=True)),
                ('created_at', models.DateTimeField(auto_now_add=True, verbose_name='creation date')),
                ('updated_at', models.DateTimeField(auto_now=True, verbose_name='last update')),
            ],
            options={
                'verbose_name': 'Category',
                'verbose_name_plural': 'Categories',
                'ordering': ['id'],
            },
        ),
        migrations.CreateModel(
            name='GeneticElement',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=255, unique=True)),
                ('glyph_thumbnail', models.FileField(upload_to=app.genetic_element.models.upload_image, validators=[app.genetic_element.models.validate_svg], verbose_name='Genetic Element Icon')),
                ('role', models.URLField(help_text='Uniform Resource Identifier (URI). Found here: http://www.sequenceontology.org/browser/obob.cgi', null=True, unique=True)),
                ('visible', models.BooleanField(default=True, help_text='Determine if it is visible in the application')),
                ('default', models.BooleanField(default=False, help_text='Determine if it will be used by default')),
                ('order', models.IntegerField(unique=True)),
                ('created_at', models.DateTimeField(auto_now_add=True, verbose_name='creation date')),
                ('updated_at', models.DateTimeField(auto_now=True, verbose_name='last update')),
                ('category', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='elements', to='genetic_element.Category')),
            ],
            options={
                'verbose_name': 'Genetic Element',
                'verbose_name_plural': 'Genetic Elements',
                'ordering': ['order'],
            },
        ),
    ]
