import os

from django.core.management.base import BaseCommand


class Command(BaseCommand):

    def handle(self, *args, **options):
        self.stdout.write('Removing migrations files...')
        self.stdout.write('Removing .py files...')
        os.system('find . -path "*/migrations/*.py" -not -name "__init__.py" -delete')
        self.stdout.write('Removing .pyc files...')
        os.system('find . -path "*/migrations/*.pyc"  -delete')
