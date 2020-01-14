import { Component, OnInit } from '@angular/core';
import { ActivatedRoute, Router } from '@angular/router';
import { finalize } from 'rxjs/operators';

import { UserHistory } from '../shared/user-history';
import { HistoryService } from '../shared/history.service';
import { JobService } from '../shared/job.service';
import { NotifyService } from 'src/app/shared/services/notify.service';
import { Job } from '../shared/job';
import { Construct } from 'src/app/construct/shared/construct';

@Component({
  selector: 'sqy-history',
  templateUrl: './history.component.html',
  styleUrls: ['./history.component.scss']
})
export class HistoryComponent implements OnInit {

  history: UserHistory = null;
  isLoading = false;
  isJobLoading = false;
  isSearching = false;
  msg: string = null;
  seq: string = null;
  resultData: { construct: Construct; results: any };

  constructor(
    private route: ActivatedRoute,
    private router: Router,
    private historySrvc: HistoryService,
    private jobSrvc: JobService,
    private notify: NotifyService
  ) { }

  ngOnInit() {
    this.history = new UserHistory().deserialize(this.route.snapshot.data.history);
    this.getJob();
  }

  getJob() {
    this.isJobLoading = true;
    this.jobSrvc.getById(this.history.job_id)
      .pipe(finalize(() => this.isJobLoading = false))
      .subscribe(
        (data: Job) => {
          this.history.job = data;
          // this.history.job.status = 'queued';
          if (this.history.isDone()) {
            this.resultData = {
              construct: this.history.construct,
              results: this.history.job.result
            };
          }
        },
        err => {
          console.log(err);
        }
      );
  }

  getStatusClass() {
    return {
      success: this.history.isDone(),
      active: this.history.isActive(),
      error: this.history.hasFailed()
    };
  }

  deleteHistory() {
    if (confirm('Are you sure?')) {
      this.isLoading = true;
      this.historySrvc.delete(this.history.id)
        .pipe(finalize(() => this.isLoading = false))
        .subscribe(() => {
          this.notify.success('History deleted!', 'bottom-right', true);
          this.router.navigate(['/workspace']);
        },
          err => this.notify.error(err.msg || 'Unable to delete history', 'bottom-right')
        );
    }
  }

}
