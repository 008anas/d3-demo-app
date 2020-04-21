import { TestBed } from '@angular/core/testing';

import { SentryErrorHandler } from './sentry-error-handler.service';

describe('SentryErrorHandlerService', () => {
  beforeEach(() => TestBed.configureTestingModule({}));

  it('should be created', () => {
    const service: SentryErrorHandler = TestBed.inject(SentryErrorHandler);
    expect(service).toBeTruthy();
  });
});
