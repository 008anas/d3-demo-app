import { Component } from '@angular/core';
import { finalize } from 'rxjs/operators';

import { NzMessageService } from 'ng-zorro-antd/message';

import { main } from '@config/main';
import { Contact } from '@models/contact';
import { SqrutinyService } from '@services/sqrutiny.service';

@Component({
  selector: 'sqy-contact',
  templateUrl: './contact.component.html',
  styleUrls: ['./contact.component.scss']
})
export class ContactComponent {

  isLoading = false;
  contact: Contact = new Contact();
  email = '';

  constructor(
    private sqrutinySrvc: SqrutinyService,
    private notify: NzMessageService
  ) {
    this.email = main.email;
  }

  submitForm() {
    this.isLoading = true;
    this.sqrutinySrvc.contact(this.contact)
      .pipe(finalize(() => this.isLoading = false))
      .subscribe(
        () => {
          // console.log('POST call successful value returned in body');
        },
        err => this.notify.error(err));
  }
}
