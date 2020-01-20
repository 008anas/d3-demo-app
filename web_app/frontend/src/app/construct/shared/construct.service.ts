import { Injectable } from '@angular/core';
import { HttpClient } from '@angular/common/http';
import { Observable } from 'rxjs';

import { environment as env } from 'src/environments/environment';
import { Construct } from './construct';

const URL_ENV = '/constructs/';

@Injectable({ providedIn: 'root' })
export class ConstructService {

  private url: string;  // URL to web api

  constructor(private http: HttpClient) {
    this.url = `${env.endpoints.api}${URL_ENV}`;
  }

  /** Gets a construct by id provided */
  getAll(): Observable<Construct[]> {
    return this.http.get<Construct[]>(`${this.url}`).pipe();
  }

  getById(id: string): Observable<Construct> {
    return this.http.get<Construct>(`${this.url}${id}`).pipe();
  }

  getExample(): Observable<Construct[]> {
    return this.http.get<Construct[]>(`${this.url}example`).pipe();
  }

}
