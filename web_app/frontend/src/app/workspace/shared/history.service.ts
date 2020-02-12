import { Injectable } from '@angular/core';
import { HttpClient } from '@angular/common/http';
import { Observable } from 'rxjs';

import { environment as env } from '@env/environment';
import { UserHistory } from './user-history';

const URL_ENV = '/workspace';

@Injectable({ providedIn: 'root' })
export class HistoryService {

  private url: string;  // URL to web api

  constructor(private http: HttpClient) {
    this.url = env.endpoints.api + URL_ENV;
  }

  getAll(): Observable<any> {
    return this.http.get<any>(`${this.url}/`).pipe();
  }

  getById(id: string): Observable<UserHistory> {
    return this.http.get<UserHistory>(`${this.url}/${id}`).pipe();
  }

  getCount(): Observable<any> {
    return this.http.get<any>(`${this.url}/count`).pipe();
  }

  delete(id: string): Observable<any> {
    return this.http.delete(`${this.url}/${id}`).pipe();
  }

  deleteAll(): Observable<any> {
    return this.http.delete(`${this.url}/all`).pipe();
  }

}
