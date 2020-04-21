import { Injectable } from '@angular/core';
import { HttpClient } from '@angular/common/http';
import { Observable } from 'rxjs';

import { environment as env } from '@env/environment';
import { Track } from './track';

const URL_ENV = '/tracks/';

@Injectable({ providedIn: 'root' })
export class TrackService {

  private url: string;  // URL to web api

  constructor(private http: HttpClient) {
    this.url = env.endpoints.api + URL_ENV;
  }

  /** Gets a track by id provided */
  getAll(): Observable<Track[]> {
    return this.http.get<Track[]>(`${this.url}`).pipe();
  }

  getByCategories(): Observable<any[]> {
    return this.http.get<any[]>(`${this.url}categories`).pipe();
  }

}
