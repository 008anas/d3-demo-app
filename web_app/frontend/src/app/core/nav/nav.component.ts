import { Component, OnInit } from '@angular/core';
import { finalize } from 'rxjs/operators';

import { HistoryService } from 'src/app/workspace/shared/history.service';

@Component({
  selector: 'sqy-main-nav',
  templateUrl: './nav.component.html',
  styleUrls: ['./nav.component.scss']
})
export class NavComponent implements OnInit {

  isLoading: boolean = false;
  count: number = 0;

  constructor(
    private historySrvc: HistoryService
  ) {
  }

  ngOnInit() {
    this.getHistoriesCount();
  }

  getHistoriesCount() {
    this.historySrvc.getCount()
      .pipe(finalize(() => this.isLoading = false))
      .subscribe(
        data => this.count = data.count,
        err => {
          console.log(err)
        }
      );
  }

}
