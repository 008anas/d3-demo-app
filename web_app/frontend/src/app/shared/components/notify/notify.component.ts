import { Component, OnInit, OnDestroy } from '@angular/core';
import { Subscription } from 'rxjs/internal/Subscription';

import { NotifyService } from '../../services/notify.service';
import { Notify, NotifyType } from '../../models/notify';

@Component({
  selector: 'notify',
  template: `<div *ngFor="let notify of notifies" class="notify {{notify.position}} do-show" [attr.data-notification-status]="cssClass(notify)">
    {{notify.message}}
    <i class="close link icon" (click)="remove(notify)"></i>
  </div>`,
  styleUrls: ['./notify.component.scss']
})
export class NotifyComponent implements OnInit, OnDestroy {

  private subscription: Subscription;
  notifies: Notify[] = [];

  constructor(private notifySrvc: NotifyService) { }

  ngOnInit() {
    this.subscription = this.notifySrvc.getMessage().subscribe(message => {
      if (!message) {
        // clear notifies when an empty notify is received
        this.notifies = [];
        return;
      }
      // add notify to array
      this.notifies.push(message);
    });
  }

  ngOnDestroy() {
    this.subscription.unsubscribe();
  }

  remove(notify: Notify) {
    this.notifies = this.notifies.filter(x => x !== notify);
  }

  cssClass(notify: Notify) {
    if (!notify) {
      return null;
    }

    // return css class based on notify type
    switch (notify.type) {
      case NotifyType.Success:
        return 'success';
      case NotifyType.Error:
        return 'error';
      case NotifyType.Info:
        return 'notice';
      case NotifyType.Warning:
        return 'warning';
    }
  }

}
