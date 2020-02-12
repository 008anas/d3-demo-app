import { Deserializable } from 'app/shared/models/deserializable.model';

export class Job extends Deserializable {
  uuid: string;
  status: string;
  result: any;
}
