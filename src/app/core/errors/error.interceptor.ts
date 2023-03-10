import { Injectable } from '@angular/core';
import { HttpRequest, HttpHandler, HttpEvent, HttpInterceptor } from '@angular/common/http';
import { catchError } from 'rxjs/operators';
import { Router } from '@angular/router';
import { Observable } from 'rxjs/internal/Observable';

import { routes } from '@config/routes';
import { environment as env } from '@env/environment';
import { throwError } from 'rxjs/internal/observable/throwError';

@Injectable()
export class ErrorInterceptor implements HttpInterceptor {

  constructor(
    private router: Router
  ) { }

  intercept(request: HttpRequest<any>, next: HttpHandler): Observable<HttpEvent<any>> {
    return next.handle(request).pipe(catchError(err => {
      if (env.production) {
        this.processError(err.status);
      }
      return throwError(err.error || err.message);
    }));
  }

  private processError(error: any) {
    switch (error) {
      case 0:
      case 500:
        this.router.navigate([routes.error500]);
        break;
      case 404:
        this.router.navigate([routes.error404]);
        break;
      default:
    }
  }
}
