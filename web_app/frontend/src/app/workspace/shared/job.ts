import { Deserializable } from '@models/deserializable.model';

export class Job extends Deserializable {
  uuid: string;
  status: string;
  result: any;
}
