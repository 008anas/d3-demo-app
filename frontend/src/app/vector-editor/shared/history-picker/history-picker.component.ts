import { Component, OnInit } from '@angular/core';
import { Subscription } from 'rxjs/internal/Subscription';
import { finalize } from 'rxjs/operators';

import { HistoryService } from 'app/workspace/shared/history.service';
import { UserHistory } from 'app/workspace/shared/user-history';

@Component({
  selector: 'sqy-history-picker',
  templateUrl: './history-picker.component.html',
  styleUrls: ['./history-picker.component.scss']
})
export class HistoryPickerComponent implements OnInit {

  sub: Subscription;
  histories: UserHistory[] = [];
  editor: any = null;
  response: any = null;
  isLoading: boolean;
  selected: UserHistory = null;

  constructor(private historySrvc: HistoryService) { }

  ngOnInit(): void {
    this.isLoading = true;
    this.getUserHistory();
  }

  getUserHistory() {
    this.isLoading = true;
    this.historySrvc.getAll()
      .pipe(finalize(() => this.isLoading = false))
      .subscribe(data => this.histories = data.map((h: any) => new UserHistory().deserialize(h)) || []);
  }

  handleSelect(history: UserHistory, state: boolean) {
    if (state) {
      this.selected = history;
    }
  }

}
