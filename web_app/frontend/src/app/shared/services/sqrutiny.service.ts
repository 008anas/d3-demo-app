import { Injectable } from '@angular/core';
import { HttpClient } from '@angular/common/http';
import { Observable } from 'rxjs';

import { environment as env } from '@env/environment';
import { Construct } from '@models/construct';
import { UserHistory } from 'app/workspace/shared/user-history';

@Injectable({ providedIn: 'root' })
export class SqrutinyService {

  private url: string;  // URL to web api

  constructor(private http: HttpClient) {
    this.url = env.endpoints.api;
  }

  fromSketch(construct: Construct): Observable<UserHistory> {
    return this.http.post<UserHistory>(`${this.url}/optimize_seq/from-sketch`, construct).pipe();
  }

  motifInSeq(sequence: string, motif: string): Observable<any> {
    if (sequence && motif) {
      return this.http.get(`${this.url}/search-motif`, { params: { sequence, motif } }).pipe();
    }
  }

  contact(contact: any) {
    return this.http.post<any[]>(`${env.endpoints.api}/contact/`, contact).pipe();
  }

}
