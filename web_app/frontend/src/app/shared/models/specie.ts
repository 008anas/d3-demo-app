import { Deserializable } from './deserializable.model';

export class Specie implements Deserializable {
  name: string;
  slug: string;
  default: boolean;
  ncbi_tax_id: number;

  deserialize(input: any) {
    Object.assign(this, input);
    return this;
  }
}
