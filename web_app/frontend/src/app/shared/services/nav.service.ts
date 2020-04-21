import { Injectable } from '@angular/core';
import { Subject } from 'rxjs/internal/Subject';
import { Observable } from 'rxjs/internal/Observable';

@Injectable({
  providedIn: 'root'
})
export class NavService {

  _badgeStatus: Subject<void> = new Subject();

  updateBadge() {
    this._badgeStatus.next();
  }

  badgeStatus(): Observable<void> {
    return this._badgeStatus.asObservable();
  }
}
