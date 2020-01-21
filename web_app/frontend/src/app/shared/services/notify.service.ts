import { Injectable } from '@angular/core';
import { Router, NavigationStart } from '@angular/router';
import { Observable, Subject } from 'rxjs';

import { NotifyType, Notify } from '../models/notify';

@Injectable({
  providedIn: 'root'
})
export class NotifyService {

  private subject = new Subject<any>();
  private keepAfterNavigation = false;

  constructor(private router: Router) {
    // clear alert message on route change
    this.router.events.subscribe(event => {
      if (event instanceof NavigationStart) {
        this.keepAfterNavigation ? this.keepAfterNavigation = false : this.subject.next();
      }
    });
  }

  success(message: string, pos = 'bottom-right', keepAfterNavigation = false) {
    this.notify(NotifyType.Success, message, pos, keepAfterNavigation);
  }

  error(message: string, pos = 'bottom-right', keepAfterNavigation = false) {
    this.notify(NotifyType.Error, message, pos, keepAfterNavigation);
  }

  info(message: string, pos = 'bottom-right', keepAfterNavigation = false) {
    this.notify(NotifyType.Info, message, pos, keepAfterNavigation);
  }

  warn(message: string, pos = 'bottom-right', keepAfterNavigation = false) {
    this.notify(NotifyType.Warning, message, pos, keepAfterNavigation);
  }

  notify(type: NotifyType, message: string, pos = 'bottom-right', keepAfterNavigation = false) {
    this.keepAfterNavigation = keepAfterNavigation;
    this.subject.next({ type: type, message: message, position: pos } as Notify);
  }

  getMessage(): Observable<any> {
    return this.subject.asObservable();
  }
}
