import { Component, OnInit, OnDestroy } from '@angular/core';
import { Subscription } from 'rxjs/internal/Subscription';

import { Notify, NotifyType } from '@models/notify';
import { NotifyService } from '@services/notify.service';

@Component({
  selector: 'sqy-notify',
  templateUrl: './notify.component.html',
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