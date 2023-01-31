import { Component } from '@angular/core';

@Component({
  selector: 'sqy-error404',
  template: `<div class="vh-100">
    <div class="centered-content">
      <h1 class="main-400">404</h1>
      <h2 class="second-h">Oops Page Not Found</h2>
      <p>The page you are looking for does not exist or has been moved.</p>
      <a href="/" class="ui primary button">Go to Home</a>
    </div>
  </div>`,
  styleUrls: ['./error404.component.scss']
})
export class Error404Component {
}
