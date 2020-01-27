# Generated by Django 2.2.8 on 2020-01-27 12:03

import uuid

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    initial = True

    dependencies = [
        ('construct', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='History',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('uuid', models.UUIDField(default=uuid.uuid4, editable=False, unique=True)),
                ('name', models.CharField(blank=True, max_length=255)),
                ('job_id', models.UUIDField()),
                ('request_ip', models.GenericIPAddressField(null=True)),
                ('deleted', models.BooleanField(default=False)),
                ('created_at', models.DateTimeField(auto_now_add=True, verbose_name='creation date')),
                ('updated_at', models.DateTimeField(auto_now=True, verbose_name='last update')),
                ('construct', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='construct.Construct')),
            ],
            options={
                'verbose_name_plural': 'Histories',
            },
        ),
    ]
