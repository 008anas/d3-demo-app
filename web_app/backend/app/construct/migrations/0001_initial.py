# Generated by Django 2.2.8 on 2020-02-10 09:52

import uuid

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    initial = True

    dependencies = [
        ('specie', '0001_initial'),
        ('genetic_element', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='Construct',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('uuid', models.UUIDField(default=uuid.uuid4, editable=False, unique=True)),
                ('name', models.CharField(max_length=255)),
                ('dna_seq', models.TextField(blank=True)),
                ('protein_seq', models.TextField(blank=True)),
                ('circular', models.BooleanField(blank=True, default=False)),
                ('description', models.TextField(blank=True, null=True)),
                ('example', models.BooleanField(default=False, help_text='Determine if is it used as example')),
                ('from_file', models.BooleanField(default=False, help_text='Describe if construct was loaded from file or using SQRUTINY sketcher')),
                ('deleted', models.BooleanField(default=False, help_text='Determine if it is visible in the application')),
                ('created_at', models.DateTimeField(auto_now_add=True, verbose_name='creation date')),
                ('updated_at', models.DateTimeField(auto_now=True, verbose_name='last update')),
                ('specie', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='specie.Specie')),
            ],
        ),
        migrations.CreateModel(
            name='Track',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('label', models.CharField(blank=True, max_length=255, null=True)),
                ('sequence', models.TextField()),
                ('start', models.IntegerField(blank=True)),
                ('end', models.IntegerField(blank=True)),
                ('color', models.CharField(default='#4e0a77', max_length=10)),
                ('construct', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='construct.Construct')),
                ('genetic_element', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='genetic_element.GeneticElement')),
            ],
        ),
        migrations.AddField(
            model_name='construct',
            name='tracks',
            field=models.ManyToManyField(through='construct.Track', to='genetic_element.GeneticElement'),
        ),
    ]
