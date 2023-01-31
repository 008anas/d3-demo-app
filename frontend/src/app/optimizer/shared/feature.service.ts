import { Injectable } from '@angular/core';
import { HttpClient } from '@angular/common/http';

import { environment as env } from '@env/environment';
import { Feature } from './feature';
import { Observable } from 'rxjs/internal/Observable';

const URL_ENV = '/features/';

@Injectable({
  providedIn: 'root'
})
export class FeatureService {

  private url: string;  // URL to web api

  constructor(private http: HttpClient) {
    this.url = env.endpoints.api + URL_ENV;
  }

  getAll(specie_id?: any): Observable<Feature[]> {
    if (specie_id) {
      return this.http.get<Feature[]>(`${this.url}`, { params: { specie_id } }).pipe();
    }
    return this.http.get<Feature[]>(`${this.url}`).pipe();
  }

}
