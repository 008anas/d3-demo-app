# Generated by Django 2.2.8 on 2020-01-29 11:33

from django.db import migrations, models


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Contact',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(blank=True, help_text='Sender name', max_length=255)),
                ('email', models.EmailField(help_text='Sender email', max_length=254)),
                ('subject', models.CharField(blank=True, help_text='Email subject', max_length=255)),
                ('message', models.TextField(help_text='Email body')),
                ('created_at', models.DateTimeField(auto_now_add=True, verbose_name='creation date')),
                ('updated_at', models.DateTimeField(auto_now=True, verbose_name='last update')),
            ],
        ),
    ]
