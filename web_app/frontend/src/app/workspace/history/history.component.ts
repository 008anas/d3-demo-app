import { Component, OnInit, OnDestroy } from '@angular/core';
import { ActivatedRoute, Router } from '@angular/router';
import { finalize } from 'rxjs/operators';

import { NzModalService } from 'ng-zorro-antd/modal';
import { NzMessageService } from 'ng-zorro-antd/message';

import { UserHistory } from '../shared/user-history';
import { HistoryService } from '../shared/history.service';
import { JobService } from '../shared/job.service';
import { Job } from '../shared/job';
import { Construct } from '@models/construct';
import { TextModalComponent } from '@components/text-modal/text-modal.component';
import { TitleService } from '@services/title.service';
import { Track } from 'app/optimizer/shared/track';

const MAX_ATTEMPTS = 15;
const RETRY_IN = 5000;

@Component({
  selector: 'sqy-history',
  templateUrl: './history.component.html',
  styleUrls: ['./history.component.scss']
})
export class HistoryComponent implements OnInit, OnDestroy {

  history: UserHistory = null;
  isLoading = false;
  isJobLoading = false;
  isSearching = false;
  msg: string = null;
  seq: string = null;
  resultData: { construct: Construct; results: any };
  trackHovered: Track = null;
  interval: any;
  attempts = 0;

  constructor(
    private route: ActivatedRoute,
    private router: Router,
    private historySrvc: HistoryService,
    private jobSrvc: JobService,
    private modal: NzModalService,
    private notify: NzMessageService,
    private titleSrvc: TitleService
  ) { }

  ngOnInit() {
    this.isLoading = true;
    this.history = new UserHistory().deserialize(this.route.snapshot.data.history);
    this.titleSrvc.setTitle(this.history.name);
    this.getJob();
  }

  ngOnDestroy() {
    if (this.interval) {
      clearInterval(this.interval);
    }
  }

  getJob() {
    this.isJobLoading = true;
    this.jobSrvc.getById(this.history.job_id)
      .pipe(finalize(() => {
        this.isLoading = false;
        this.isJobLoading = false;
      }))
      .subscribe(
        (data: Job) => {
          this.history.job = data;
          // this.history.job.status = 'queued';
          if (this.history.isDone()) {
            if (this.interval) {
              clearInterval(this.interval);
            }
            this.resultData = {
              construct: this.history.construct,
              results: this.history.job.result
            };
          } else if (this.history.isActive()) {
            if (!this.interval) {
              this.interval = setInterval(() => this.getJob(), RETRY_IN);
              this.attempts = 0;
            }
            this.attempts++;
            if (this.attempts >= MAX_ATTEMPTS) { clearInterval(this.interval); }
          }
        },
        err => {
          console.log(err);
        }
      );
  }

  deleteHistory() {
    if (confirm('Are you sure?')) {
      this.isLoading = true;
      this.historySrvc.delete(this.history.id)
        .pipe(finalize(() => this.isLoading = false))
        .subscribe(() => {
          this.notify.success('History deleted!');
          this.router.navigate(['/workspace']);
        },
          err => this.notify.error(err || 'Unable to delete history')
        );
    }
  }

  textModal(str: string) {
    this.modal.create({
      nzContent: TextModalComponent,
      nzWrapClassName: 'center-modal',
      nzComponentParams: {
        txt: str
      },
      nzFooter: null
    });
  }

}
