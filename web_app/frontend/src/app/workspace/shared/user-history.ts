import { Deserializable } from 'app/shared/models/deserializable.model';
import { Job } from './job';
import { Construct } from '@models/construct';

export class UserHistory extends Deserializable {
  id: string;
  name: string;
  slug: string;
  construct: Construct;
  job_id: string;
  job: Job;
  created_at: string;
  updated_at: string;

  isDone() {
    return this.job && this.job.status === 'finished';
  }

  isActive() {
    return this.job && (this.job.status === 'queued' || this.job.status === 'started');
  }

  hasFailed() {
    return !this.isActive() && !this.isDone();
  }
}
