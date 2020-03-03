import { Component, OnInit } from '@angular/core';
import { finalize } from 'rxjs/operators';

import { main } from '@config/main';
import { Contact } from '@models/contact';
import { SqrutinyService } from '@services/sqrutiny.service';

@Component({
  selector: 'sqy-contact',
  templateUrl: './contact.component.html',
  styleUrls: ['./contact.component.scss']
})
export class ContactComponent implements OnInit {

  isLoading = false;
  contact: Contact = new Contact();
  email: string;

  constructor(private sqrutinySrvc: SqrutinyService) { }

  ngOnInit() {
    this.email = main.email;
  }

  submitForm() {
    this.isLoading = true;
    this.sqrutinySrvc.contact(this.contact)
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
