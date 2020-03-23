import { Deserializable } from '@models/deserializable.model';
import { Track } from 'app/optimizer/shared/track';
import { Specie } from '@models/specie';

export class Construct extends Deserializable {
  id?: string;
  name: string;
  circular?: boolean;
  description: string;
  dna_seq: string;
  protein_seq: string;
  specie: Specie;
  specie_tax_id: number;
  tracks: Track[];
  n_tracks: number;
}
