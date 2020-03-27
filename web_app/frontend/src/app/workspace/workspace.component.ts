import { Component, OnInit } from '@angular/core';
import { Router } from '@angular/router';
import { finalize } from 'rxjs/operators';

import { NzMessageService } from 'ng-zorro-antd/message';

import { HistoryService } from './shared/history.service';
import { UserHistory } from './shared/user-history';
import { NavService } from '@services/nav.service';

@Component({
  selector: 'sqy-workspace',
  templateUrl: './workspace.component.html',
  styleUrls: ['./workspace.component.scss']
})
export class WorkspaceComponent implements OnInit {

  histories: UserHistory[] = [];
  isLoading = false;

  constructor(
    private router: Router,
    private historySrvc: HistoryService,
    private notify: NzMessageService,
    private navSrvc: NavService
  ) {
    this.router.routeReuseStrategy.shouldReuseRoute = () => false;
  }

  ngOnInit() {
    this.isLoading = true;
    this.getUserHistory();
  }

  getUserHistory() {
    this.isLoading = true;
    this.historySrvc.getAll()
      .pipe(finalize(() => this.isLoading = false))
      .subscribe(data => this.histories = data.map((h: any) => new UserHistory().deserialize(h)) || []);
  }

  clearHistory() {
    if (confirm('Are you sure that you want to clear history?. This action cannot be reverted.')) {
      this.isLoading = true;
      this.historySrvc.deleteAll()
        .pipe(finalize(() => this.isLoading = false))
        .subscribe(data => {
          this.navSrvc.updateBadge();
          this.getUserHistory();
          this.notify.success(data.msg);
        });
    }
  }

  deleteHistory(id: string) {
    if (confirm('Are you sure that you want to delete such history?. This action cannot be reverted.')) {
      this.isLoading = true;
      this.historySrvc.delete(id)
        .pipe(finalize(() => this.isLoading = false))
        .subscribe(() => {
          this.navSrvc.updateBadge();
          this.getUserHistory();
          this.notify.success('History deleted!');
        });
    }
  }

}
