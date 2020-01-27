import { Injectable } from '@angular/core';
import { HttpRequest, HttpHandler, HttpEvent, HttpInterceptor } from '@angular/common/http';
import { Observable, throwError } from 'rxjs';
import { catchError } from 'rxjs/operators';
import { Router } from '@angular/router';

import { environment as env } from 'src/environments/environment';
import { routes } from 'src/app/config/routes';

@Injectable()
export class ErrorInterceptor implements HttpInterceptor {

  constructor(
    private router: Router
  ) { }

  intercept(request: HttpRequest<any>, next: HttpHandler): Observable<HttpEvent<any>> {
    return next.handle(request).pipe(catchError(err => {
      this.processError(err.status);
      return throwError(err.error || err.statusText);
    }));
  }

  private processError(error: any) {
    if (env.production) {
      switch (error) {
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
}
