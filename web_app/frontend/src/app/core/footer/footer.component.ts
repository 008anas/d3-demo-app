import { Component } from '@angular/core';

@Component({
  selector: 'sqy-footer',
  templateUrl: './footer.component.html',
  styleUrls: ['./footer.component.scss']
})
export class FooterComponent {

  public getCurrentYear(){
    return new Date().getFullYear();
  }
}
