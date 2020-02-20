import { Component } from '@angular/core';

@Component({
  selector: 'sqy-error500',
  template: `<div class="vh-100">
    <div class="centered-content">
      <h1 class="main-500">500</h1>
      <h2 class="second-h">Internal Server Error</h2>
      <p>Something goes wrong in the server, please contact <strong>dbspipes@crg.es</strong> to get support or try again later.</p>
      <a href="/" class="ui primary button">Return to Homepage</a>
    </div>
  </div>`,
  styleUrls: ['./error500.component.scss']
})
export class Error500Component {
}
