import { Deserializable } from '@models/deserializable.model';

export class Feature extends Deserializable {
  id?: number;
  name: string;
  alias: string;
  description: string;
  genome_min: number;
  genome_max: number;
}
