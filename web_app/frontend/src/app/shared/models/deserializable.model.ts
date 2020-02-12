import { DeserializableI } from 'app/shared/models/deserializable.interface';

export class Deserializable implements DeserializableI {

  deserialize(input: any) {
    Object.assign(this, input);
    return this;
  }
}
