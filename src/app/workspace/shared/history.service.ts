import { Injectable } from '@angular/core';
import { HttpClient } from '@angular/common/http';
import { Observable } from 'rxjs/internal/Observable';

import { environment as env } from '@env/environment';
import { UserHistory } from './user-history';
import { of } from 'rxjs';

const URL_ENV = '/workspace';

@Injectable({ providedIn: 'root' })
export class HistoryService {

  private url: string;  // URL to web api

  constructor(private http: HttpClient) {
    this.url = env.endpoints.api + URL_ENV;
  }

  getAll(): Observable<any> {
  return of([ {id: 'ecbc2279-0d4e-4844-83b4-e876ef089929',
    name: 'Example Job',
    slug: 'example-job',
    construct: {},
    job_id: 1,
    status: 'finished',
    job: {
      status: 'finished'
    },
    created_at: '01/01/2023',
    updated_at: '0101010'}]);

    // return this.http.get<any>(`${this.url}/`).pipe();
  }

  getById(id: string): Observable<any> {
    return of( {id: 'ecbc2279-0d4e-4844-83b4-e876ef089929',
      name: 'Example Job',
      slug: 'example-job',
      construct: {
        dna_seq: 'GACTCTGCCGCAGCCACGTATCGCCTGAAAGCCAGTCAGCGTTAAGGAGTGCTCTGAGCAGGACAACTCGCGTAGTGAGAGTTACATGTTCGTTGGGCTCTTCCGACACGGACCTGAGTTGGCCAACGTCCCACCTGAGGTCTGTGCCCCGGTGATGAGAAGTGTGCATCTCGTTCTTGCAGCTCGTCAGTACTTTCAGAATCATGGCGTGCATGGTAGAATGACCCTTATAACGGACTTCGACATGGCAATAACCCCCCGTTTCTACTTCTAGAGGAGAAAAGTATTGACATGAGCGCTCCCGGCACAAGGGCCAAAGAAGTCTCCAATTTCTTATTTCCGAATGACATGCGTCTCCTTGCGGGTAAATCACCGACCGCAATTCATAGAAGCCTGGGGGAACAGATAGGTCTAATTAGCTTAAGAGAGTAAATCCTGGGATCATTCAGTAGTAACCATAAACTTACGCTGGGGCTTCTTCGGCGGATTTTTACAGTTAC',
        protein_seq: 'IVYFARMKHSDDSWMQYHTRWTWACLTHAATCSAMYLMADIEPHFWSEHPQCIQTRYHPVVPETHESNKPVTANNLYKANNMSHPPDAYAQIKWCFQGNNARFDAPFEGRWGHIPMVEITGTLHSYAQNWIIHGWIVIRHNGNPCMTDEYMNTWGKEIIDDDINHNLLGEFEYMVLFVPWYTFLMPVQETMFNHARENRDNWHVFRIVSCTMYATQRIPECVQWMTSDYDMNCFSCSTQYSHMCMLLKLLGRIKCDFAWRRRYCWWVLPHDYHMVLMINQHQDLYQKRQNGPGNVCGCDVRMKVRVLPELPPNESKTNLCYWARFMNYRKSKAWEYCSWTDPYESAEWYFFGWYAWPVVFLTVQPPFAVGPMKNGFMREMDYIRAASQHAEVRMPNNSPELINYNVIRPYKQKEMYQAYYLIILVYWVWSNDVQHWGTIQHAWIIDLIKPCRCRKAEFHAVTLYKRCFDLKETKTLFESYWHGSQTQDWGEPWTLKWVIW',
        tracks: [{
          id: 1,
          type: 'string',
          label: 'string',
          color: 'red',
          start: 0,
          end: 50,
          pos: 1,
          sequence: 'GACTCTG'
        }, {
          id: 1,
          type: 'string',
          label: 'string',
          color: 'yellow',
          start: 50,
          end: 190,
          pos: 2,
          sequence: 'GACTCTG'
        }, {
          id: 1,
          type: 'string',
          label: 'string',
          color: 'blue',
          start: 314,
          end: 459,
          pos: 3,
          sequence: 'GACTCTG'
        }]
      },
      job_id: 1,
      status: 'finished',
      job: {
        status: 'finished'
      },
      created_at: '01/01/2023',
      updated_at: '0101010'});
    // return this.http.get<UserHistory>(`${this.url}/${id}`).pipe();
  }

  getByIdNot404(id: string): Observable<UserHistory> {
    return this.http.get<UserHistory>(`${this.url}/history/${id}`).pipe();
  }

  getCount(): Observable<any> {
    return of({count: 1});
  }

  delete(id: string): Observable<any> {
    return this.http.delete(`${this.url}/${id}`).pipe();
  }

  deleteAll(): Observable<any> {
    return this.http.delete(`${this.url}/all`).pipe();
  }

  export(id: string, params?) {
    return this.http.post(`${this.url}/export/${id}`, params).pipe();
  }

}
