import { Injectable } from '@angular/core';
import { HttpClient } from '@angular/common/http';

import { environment as env } from '@env/environment';
import { Feature } from './feature';
import { Observable } from 'rxjs/internal/Observable';
import { of } from 'rxjs';

const URL_ENV = '/features/';

@Injectable({
  providedIn: 'root'
})
export class FeatureService {

  private url: string;  // URL to web api

  constructor(private http: HttpClient) {
    this.url = env.endpoints.api + URL_ENV;
  }

  getAll(specie_id?: any): Observable<any[]> {
    return of([
      { id: 1,
        name: 'Specie 1',
        alias: 'Specie 1',
        description: 'Test',
        genome_min: 1,
        genome_max: 10,}
    ])
    // if (specie_id) {
    //   return this.http.get<Feature[]>(`${this.url}`, { params: { specie_id } }).pipe();
    // }
    // return this.http.get<Feature[]>(`${this.url}`).pipe();
  }

}
