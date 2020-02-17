import { Injectable } from '@angular/core';
import { ActivatedRouteSnapshot, Resolve } from '@angular/router';
import { Observable } from 'rxjs';

import { Construct } from '@models/construct';
import { ConstructService } from '@services/construct.service';

@Injectable()
export class ConstructResolver implements Resolve<Observable<Construct>> {

  constructor(private constructSrvc: ConstructService) {
  }

  resolve(route: ActivatedRouteSnapshot) {
    return this.constructSrvc.getById(route.paramMap.get('construct'));
  }
}
