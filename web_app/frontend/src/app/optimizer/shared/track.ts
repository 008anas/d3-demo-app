import { Deserializable } from 'src/app/shared/models/deserializable.model';

export class Track extends Deserializable {
  id?: number;
  name: string;
  label?: string;
  color: string;
  start?: string;
  end?: string;
  sequence: string;
  default?: boolean;
  selected?: boolean;
  glyph_thumbnail: string;
}
