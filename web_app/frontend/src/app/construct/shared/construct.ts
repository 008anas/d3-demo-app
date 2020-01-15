import { Deserializable } from 'src/app/shared/models/deserializable.model';
import { Track } from 'src/app/optimizer/shared/track';
import { Specie } from 'src/app/shared/models/specie';

export class Construct extends Deserializable {
  id?: string;
  label: string;
  name?: string;
  circular?: boolean;
  sequence: string;
  dna_seq: string;
  protein_seq: string;
  specie: Specie;
  tracks: Track[];
  n_tracks: number;
}
