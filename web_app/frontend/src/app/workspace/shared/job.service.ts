import { Injectable } from '@angular/core';
import { HttpClient } from '@angular/common/http';
import { Observable } from 'rxjs';

import { environment as env } from '@env/environment';
import { Job } from './job';

const URL_ENV = '/workspace/job/';

@Injectable({ providedIn: 'root' })
export class JobService {

  private url: string;  // URL to web api

  constructor(private http: HttpClient) {
    this.url = env.endpoints.api + URL_ENV;
  }

  getById(id: string): Observable<Job> {
    return this.http.get<Job>(`${this.url}${id}`).pipe();
  }

}
