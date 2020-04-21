import { Component, OnInit, OnDestroy } from '@angular/core';
import { Subscription } from 'rxjs/internal/Subscription';

import { HistoryService } from 'app/workspace/shared/history.service';
import { NavService } from '@services/nav.service';

@Component({
  selector: 'sqy-main-nav',
  templateUrl: './nav.component.html',
  styleUrls: ['./nav.component.scss']
})
export class NavComponent implements OnInit, OnDestroy {

  count = 0;
  sub: Subscription;

  constructor(
    private historySrvc: HistoryService,
    private navSrvc: NavService
  ) { }

  ngOnInit() {
    this.getHistoriesCount();
    this.sub = this.navSrvc.badgeStatus()
      .subscribe(
        () => this.getHistoriesCount()
      );
  }

  ngOnDestroy() {
    this.sub.unsubscribe();
  }

  getHistoriesCount() {
    this.historySrvc.getCount()
      .subscribe(data => this.count = Number.parseInt(data.count, 10));
  }

}
