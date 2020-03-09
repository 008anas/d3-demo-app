import json
import os
from datetime import datetime

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
                Specie(id=s.get('id'), name=s.get('name'), comment=s.get('comment', None),
                       genome_gbk=s.get('genome_gbk', None), slug=s.get('slug'), tax_id=s.get('tax_id'),
                       tax_link=s.get('tax_link'),
                       gc_content=s.get('gc_content'), codon_table=s.get('codon_table'), default=s.get('default', None),
                       visible=s.get('visible'),
                       created_at=s.get('created_at', datetime.now()),
                       updated_at=s.get('updated_at', datetime.now())).save()
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
                Category(id=c.get('id'), name=c.get('name'), created_at=c.get('created_at', datetime.now()),
                         updated_at=c.get('updated_at', datetime.now())).save()
            except Exception as e:
                raise CommandError('Error while feeding: ' + str(e))

        self.stdout.write(
            self.style.SUCCESS('Category table was feeded successfully. ' + str(count) + ' records added.'))

    def feed_genetic_elements(self):
        count = 0
        json_data = self.open_json_file(os.path.join(self.path, '../../data/genetic_elements.json'))

        for e in json_data:
            count += 1
            try:
                GeneticElement(id=e.get('id'), role=e.get('role'), name=e.get('name'),
                               glyph_thumbnail=e.get('glyph_thumbnail'),
                               order=e.get('order'), visible=e.get('visible'),
                               default=e.get('default'),
                               category_id=e.get('category'), created_at=e.get('created_at', datetime.now()),
                               updated_at=e.get('updated_at', datetime.now())).save()
            except Exception as e:
                raise CommandError('Error while feeding: ' + str(e))

        self.stdout.write(
            self.style.SUCCESS('Genetic Elements table was feeded successfully. ' + str(count) + ' records added.'))

    def feed_parameters(self):
        count = 0
        json_data = self.open_json_file(os.path.join(self.path, '../../data/parameters.json'))
        for m in json_data:
            count += 1
            try:
                Parameter(id=m.get('id'), name=m.get('name'), alias=m.get('alias'), specie_id=m.get('specie'),
                          genetic_element_id=m.get('genetic_element'),
                          matrix_file=m.get('matrix_file'), genome_min=m.get('genome_min'),
                          genome_max=m.get('genome_max'), active=m.get('active'),
                          created_at=m.get('created_at', datetime.now()),
                          updated_at=m.get('updated_at', datetime.now())).save()
            except Exception as e:
                raise CommandError('Error while feeding: ' + str(e))

        self.stdout.write(
            self.style.SUCCESS('Parameter table was feeded successfully. ' + str(count) + ' records added.'))

    def feed(self):
        self.stdout.write('Feeding...')
        self.feed_species()
        self.feed_categories()
        self.feed_genetic_elements()
        self.feed_parameters()
        self.stdout.write(
            self.style.SUCCESS('Database was successfully feeded.'))

    def handle(self, *args, **options):
        try:
            self.feed()

        except FileNotFoundError as f:
            raise CommandError(str(f))
