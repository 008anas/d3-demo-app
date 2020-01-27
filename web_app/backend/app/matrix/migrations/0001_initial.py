# Generated by Django 2.2.8 on 2020-01-27 12:03

import django.db.models.deletion
from django.db import migrations, models

import app.matrix.models


class Migration(migrations.Migration):

    initial = True

    dependencies = [
        ('specie', '0001_initial'),
        ('genetic_element', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='Matrix',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(help_text='Used for interface display', max_length=255)),
                ('alias', models.CharField(help_text='For internal usage. Must match checker list', max_length=255)),
                ('matrix_file', models.FileField(help_text='Matrix file', upload_to=app.matrix.models.species_dir)),
                ('genome_min', models.FloatField(null=True)),
                ('genome_max', models.FloatField(null=True)),
                ('active', models.BooleanField(default=True, help_text='Determine if it can be used')),
                ('created_at', models.DateTimeField(auto_now_add=True, verbose_name='creation date')),
                ('updated_at', models.DateTimeField(auto_now=True, verbose_name='last update')),
                ('genetic_element', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='genetic_element.GeneticElement')),
                ('specie', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='specie.Specie')),
            ],
            options={
                'verbose_name': 'Matrix',
                'verbose_name_plural': 'Matrices',
                'ordering': ['created_at'],
                'unique_together': {('alias', 'specie', 'genetic_element')},
            },
        ),
    ]
