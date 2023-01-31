import { Deserializable } from './deserializable.model';

export class Specie implements Deserializable {
  name: string;
  slug: string;
  default: boolean;
  tax_id: number;
  tax_link: string;
  comment: string;

  deserialize(input: any) {
    Object.assign(this, input);
    return this;
  }
}
