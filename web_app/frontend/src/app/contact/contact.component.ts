import { Component, OnInit } from '@angular/core';
import { finalize } from 'rxjs/operators';

import { main } from '../config/main';
import { environment as env } from 'src/environments/environment';
import { Contact } from '../shared/models/contact';
import { HttpClient } from '@angular/common/http';

@Component({
  selector: 'sqy-contact',
  templateUrl: './contact.component.html',
  styleUrls: ['./contact.component.scss']
})
export class ContactComponent implements OnInit {

  isLoading = false;
  contact: Contact = new Contact();
  email: string;

  constructor(private http: HttpClient) { }

  ngOnInit() {
    this.email = main.email;
  }

  /**
  * Process the form. Send to email
  */
  submitForm() {
    this.isLoading = true;
    this.http.post<any[]>(`${env.endpoints.api}/contact/`, this.contact)
      .pipe(finalize(() => this.isLoading = false))
      .subscribe(
        () => {
          console.log('POST call successful value returned in body');
        },
        () => {
          console.log('POST call in error');
        });
  }
}
