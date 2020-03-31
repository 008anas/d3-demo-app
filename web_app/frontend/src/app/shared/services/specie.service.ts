import { Injectable } from '@angular/core';
import { HttpClient } from '@angular/common/http';
import { Observable } from 'rxjs';

import { environment as env } from '@env/environment';
import { Specie } from '@models/specie';

const URL_ENV = '/species';

@Injectable({ providedIn: 'root' })
export class SpecieService {

  private url: string;  // URL to web api

  constructor(private http: HttpClient) {
    this.url = env.endpoints.api + URL_ENV;
  }

  /** Gets a specie by id provided */
  getAll(): Observable<Specie[]> {
    return this.http.get<Specie[]>(`${this.url}/`).pipe();
  }

  /** Gets a specie by id provided */
  getBySlug(specie: Specie): Observable<Specie> {
    return this.http.get<Specie>(`${this.url}/${specie.slug}`).pipe();
  }

}
