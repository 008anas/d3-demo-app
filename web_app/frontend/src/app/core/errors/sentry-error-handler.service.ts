import { Injectable, ErrorHandler } from '@angular/core';

import * as Sentry from '@sentry/browser';

import { environment as env } from '@env/environment';

Sentry.init({
  enabled: env.production,
  dsn: env.sentry.dsn
});

@Injectable({
  providedIn: 'root'
})
export class SentryErrorHandler implements ErrorHandler {

  constructor() { }

  handleError(error: { originalError: any; }) {
    const eventId = Sentry.captureException(error.originalError || error);
    Sentry.showReportDialog({ eventId });
  }
}
