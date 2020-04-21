import json
import os

from django.core.management.base import BaseCommand, CommandError

from app.genetic_element.models import GeneticElement, Category
from app.parameter.models import Parameter
from app.specie.models import Specie


class Command(BaseCommand):
    help = 'Feeds database from json file'
    path = os.path.dirname(__file__)

    def open_json_file(self, path):
        try:
            content = open(path, 'r')
            data = json.load(content)
            return data

        except FileNotFoundError as f:
            raise CommandError(str(f))

    def feed_species(self):
        count = 0
        json_data = self.open_json_file(os.path.join(self.path, '../../data/species.json'))
        for s in json_data:
            count += 1
            try:
                Specie(**s).save()
            except Exception as e:
                raise CommandError('Error while feeding: ' + str(e))

        self.stdout.write(
            self.style.SUCCESS('Species table was feeded successfully. ' + str(count) + ' records added.'))

    def feed_categories(self):
        count = 0
        json_data = self.open_json_file(os.path.join(self.path, '../../data/categories.json'))
        for c in json_data:
            count += 1
            try:
                Category(**c).save()
            except Exception as e:
                raise CommandError('Error while feeding: ' + str(e))

        self.stdout.write(
            self.style.SUCCESS('Category table was feeded successfully. ' + str(count) + ' records added.'))

    def feed_genetic_elements(self):
        count = 0
        json_data = self.open_json_file(os.path.join(self.path, '../../data/genetic_elements.json'))

        for ge in json_data:
            count += 1
            try:
                GeneticElement(**ge).save()
            except Exception as e:
                raise CommandError('Error while feeding: ' + str(ge))

        self.stdout.write(
            self.style.SUCCESS('Genetic Elements table was feeded successfully. ' + str(count) + ' records added.'))

    def feed_parameters(self):
        count = 0
        json_data = self.open_json_file(os.path.join(self.path, '../../data/parameters.json'))
        for p in json_data:
            count += 1
            try:
                Parameter(**p).save()
            except Exception as e:
                raise CommandError('Error while feeding: ' + str(p))

        self.stdout.write(
            self.style.SUCCESS('Parameter table was feeded successfully. ' + str(count) + ' records added.'))

    def feed(self):
        self.stdout.write('Feeding...')
        self.feed_species()
        self.feed_categories()
        self.feed_genetic_elements()
        self.feed_parameters()
        self.stdout.write(
            self.style.SUCCESS('Database was successfully feeded.')
        )

    def handle(self, *args, **options):
        try:
            self.feed()

        except FileNotFoundError as f:
            raise CommandError(str(f))
