import sys

from rest_framework import serializers
from rest_framework.exceptions import ValidationError

from app.genetic_element.models import GeneticElement
from app.specie.serializers import SpecieSerializer
from sqrutiny.settings import BASE_DIR

sys.path.insert(0, BASE_DIR + '/../../dev')
from tools import is_dna_seq_valid, seq_translator

from .models import Construct, Track


class TrackRetrieveSerializer(serializers.ModelSerializer):
    glyph_thumbnail = serializers.SerializerMethodField()
    type = serializers.SerializerMethodField()

    def get_glyph_thumbnail(self, track):
        return self.context['request'].build_absolute_uri(track.genetic_element.glyph_thumbnail.url)

    def get_type(self, track):
        return track.genetic_element.name

    class Meta:
        model = Track
        exclude = ['id', 'construct', 'genetic_element']


class TrackCreateSerializer(serializers.ModelSerializer):
    type = serializers.CharField()

    def validate_sequence(self, seq):
        seq = seq.strip().replace('\n', '').upper()
        if not is_dna_seq_valid(seq):
            raise ValidationError("Sequence not valid")
        return seq

    class Meta:
        model = Track
        exclude = ['construct', 'genetic_element']


class ConstructCreateSerializer(serializers.ModelSerializer):
    specie_tax_id = serializers.IntegerField(write_only=True)
    tracks = TrackCreateSerializer(many=True)

    class Meta:
        model = Construct
        fields = ['name', 'tracks', 'description', 'specie_tax_id', 'circular']

    def create(self, validated_data):
        validated_data.pop('specie_tax_id')
        tracks = validated_data.pop('tracks')
        seq = ''.join([track['sequence'] for track in tracks])
        construct = Construct.objects.create(**validated_data, dna_seq=seq, protein_seq=seq_translator(seq))
        i = 1
        for track in tracks:
            type = track.pop('type').lower()
            element = GeneticElement.objects.filter(name__iexact=type).first()
            if element:
                if 'start' not in track or 'end' not in track:
                    track['start'] = i
                    track['end'] = i + len(track['sequence']) - 1
                    i += len(track['sequence'])
                Track.objects.create(**track, construct=construct, genetic_element=element)
            else:
                construct.delete()
                raise ValidationError("Element type: " + type + " not valid. Please enter a valid one.")
        return construct


class ConstructRetrieveSerializer(serializers.ModelSerializer):
    id = serializers.SerializerMethodField()
    specie = SpecieSerializer(read_only=True)
    tracks = serializers.SerializerMethodField()
    n_tracks = serializers.SerializerMethodField()

    def get_id(self, construct):
        return construct.uuid

    def get_tracks(self, construct):
        tracks = Track.objects.filter(construct=construct)
        track_ser = TrackRetrieveSerializer(tracks, context={'request': self.context['request']}, many=True)
        return track_ser.data

    def get_n_tracks(self, construct):
        return construct.tracks.count()

    class Meta:
        model = Construct
        exclude = ['uuid', 'example', 'from_file','deleted', 'updated_at']


class SearchMotifSerializer(serializers.Serializer):
    motif = serializers.CharField()
    sequence = serializers.CharField()

    def validate_sequence(self, value):
        value = value.upper()
        if not is_dna_seq_valid(value):
            raise ValidationError("Sequence not valid")
        return value
