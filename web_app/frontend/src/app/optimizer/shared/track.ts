import { Deserializable } from '@models/deserializable.model';

export class Track extends Deserializable {
  id?: number;
  type: string;
  label?: string;
  color: string;
  start?: number;
  end?: number;
  sequence: string;
  default?: boolean;
  selected?: boolean;
  glyph_thumbnail: string;
}
