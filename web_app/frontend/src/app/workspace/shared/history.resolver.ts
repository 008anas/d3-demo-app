import { Injectable } from '@angular/core';
import { ActivatedRouteSnapshot, Resolve } from '@angular/router';
import { Observable } from 'rxjs';

import { UserHistory } from './user-history';
import { HistoryService } from './history.service';

@Injectable()
export class HistoryResolver implements Resolve<Observable<UserHistory>> {

  constructor(private historySrvc: HistoryService) {
  }

  resolve(route: ActivatedRouteSnapshot) {
    return this.historySrvc.getById(route.paramMap.get('history'));
  }
}
