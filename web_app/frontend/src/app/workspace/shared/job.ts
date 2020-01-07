import { Deserializable } from 'src/app/shared/models/deserializable.model';

export class Job extends Deserializable {
  uuid: string;
  status: string;
  result: any;
}
