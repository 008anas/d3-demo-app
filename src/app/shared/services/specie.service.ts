import { Injectable } from '@angular/core';
import { HttpClient } from '@angular/common/http';
import { Observable } from 'rxjs/internal/Observable';

import { environment as env } from '@env/environment';
import { Specie } from '@models/specie';
import { of } from 'rxjs';

const URL_ENV = '/species/';

@Injectable({ providedIn: 'root' })
export class SpecieService {

  private url: string;  // URL to web api

  constructor(private http: HttpClient) {
    this.url = env.endpoints.api + URL_ENV;
  }

  /** Gets a specie by id provided */
  getAll(): Observable<any[]> {
    return of([
      {name: 'Mycoplasma',slug: 'mycoplasma',default: true,tax_id: 1,tax_link: '',comment: 'Test'}
    ]);
  }

  /** Gets a specie by id provided */
  getBySlug(specie: Specie): Observable<any> {
    return of(
      {name: 'Mycoplasma',slug: 'mycoplasma',default: true,tax_id: 1,tax_link: '',comment: 'Test'}
    );
    // return this.http.get<Specie>(`${this.url}${specie.slug}`).pipe();
  }

}
