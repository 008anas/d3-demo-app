# Generated by Django 2.2.8 on 2020-03-03 11:03

import app.specie.models
from django.db import migrations, models


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Specie',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=255, unique=True)),
                ('comment', models.TextField(blank=True, help_text='Short comment to display near name', null=True)),
                ('slug', models.CharField(help_text='Short name displayed in the url', max_length=255, unique=True)),
                ('tax_id', models.IntegerField(unique=True)),
                ('tax_link', models.URLField(unique=True)),
                ('gc_content', models.FloatField(blank=True, null=True)),
                ('codon_table', models.IntegerField(default=11)),
                ('genome_gbk', models.FileField(null=True, upload_to=app.specie.models.upload_genome_gbk, verbose_name='Genome GenBank file')),
                ('default', models.BooleanField(default=False, help_text='Determine which it will be used by default')),
                ('visible', models.BooleanField(default=True, help_text='Determine if it is visible in the app')),
                ('created_at', models.DateTimeField(auto_now_add=True, verbose_name='creation date')),
                ('updated_at', models.DateTimeField(auto_now=True, verbose_name='last update')),
            ],
        ),
    ]
