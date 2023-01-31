import { Component } from '@angular/core';

import { main } from '@config/main';

@Component({
  selector: 'sqy-error500',
  template: `<div class="vh-100">
    <div class="centered-content">
      <h1 class="main-500">500</h1>
      <h2 class="second-h">Internal Server Error</h2>
      <p>Something goes wrong in the server, please try again later.</p>
      <p>If problem persist please contact us thru <a routerLink="/contact">contact section</a>
      or email us to <strong>{{email}}</strong> to get support.</p>
      <a href="/" class="ui primary button">Return to Homepage</a>
    </div>
  </div>`,
  styleUrls: ['./error500.component.scss']
})
export class Error500Component {

  email = '';

  constructor() {
    this.email = main.email;
  }

}
