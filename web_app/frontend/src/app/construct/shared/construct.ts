import { Deserializable } from 'app/shared/models/deserializable.model';
import { Track } from 'app/optimizer/shared/track';
import { Specie } from 'app/shared/models/specie';

export class Construct extends Deserializable {
  id?: string;
  name: string;
  circular?: boolean;
  sequence: string;
  description: string;
  dna_seq: string;
  protein_seq: string;
  specie: Specie;
  tracks: Track[];
  n_tracks: number;
}
