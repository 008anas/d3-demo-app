import { DeserializableI } from '@models/deserializable.interface';

export class Deserializable implements DeserializableI {

  deserialize(input: any) {
    Object.assign(this, input);
    return this;
  }
}
