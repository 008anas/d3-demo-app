import { Injectable, ErrorHandler } from '@angular/core';

import * as Sentry from '@sentry/browser';

import { environment as env } from '@env/environment';

Sentry.init({
  dsn: env.sentry.dsn
});

@Injectable({
  providedIn: 'root'
})
export class SentryErrorHandler implements ErrorHandler {

  constructor() { }

  handleError(error: { originalError: any; }) {
    if (env.production) {
      Sentry.captureException(error.originalError || error);
    } else {
      console.error(error.originalError || error);
    }
  }
}
