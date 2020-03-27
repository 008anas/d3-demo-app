import { Injectable } from '@angular/core';
import { HttpRequest, HttpHandler, HttpEvent, HttpInterceptor } from '@angular/common/http';
import { Observable } from 'rxjs/internal/Observable';

import { CookieService } from 'ngx-cookie-service';

@Injectable()
export class AppHttpInterceptor implements HttpInterceptor {

  constructor(private cookieSrvc: CookieService) { }

  intercept(request: HttpRequest<any>, next: HttpHandler): Observable<HttpEvent<any>> {

    request = request.clone({
      withCredentials: true
    });

    if (!request.headers.has('Accept')) {
      request = request.clone({ headers: request.headers.set('Accept', 'application/json') });
    }

    if (this.cookieSrvc.get('csrftoken')) {
      request = request.clone({ headers: request.headers.set('X-CSRFToken', this.cookieSrvc.get('csrftoken')) });
    }

    return next.handle(request);
  }
}
