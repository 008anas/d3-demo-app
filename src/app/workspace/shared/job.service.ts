import { Construct } from './../../shared/models/construct';
import { Injectable } from '@angular/core';
import { HttpClient } from '@angular/common/http';
import { Observable } from 'rxjs/internal/Observable';

import { environment as env } from '@env/environment';
import { Job } from './job';
import { of } from 'rxjs';

const URL_ENV = '/workspace/job/';

@Injectable({ providedIn: 'root' })
export class JobService {

  private url: string;  // URL to web api

  constructor(private http: HttpClient) {
    this.url = env.endpoints.api + URL_ENV;
  }

  getById(id: string): Observable<any> {
    return of({
      uuid: 'ecbc2279-0d4e-4844-83b4-e876ef089923',
  status: 'finished',
  result: [{
    alias: 'A',
    scores: [{start: 1, raw_score: 20},
      {start: 20, raw_score: 0.80},
      {start: 20, raw_score: 0.80},
      {start: 20, raw_score: 0.80},
      {start: 20, raw_score: 0.80},
      {start: 20, raw_score: 0.80},
      {start: 20, raw_score: 0.80},
      {start: 20, raw_score: 0.80},
    ]
  }]
    })
    // return this.http.get<Job>(`${this.url}${id}`).pipe();
  }

}
