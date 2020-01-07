import { Component, Input } from '@angular/core';

@Component({
  selector: 'sqy-sidebar',
  templateUrl: './sidebar.component.html',
  styleUrls: ['./sidebar.component.scss']
})
export class SidebarComponent {

  @Input() show = false;
  @Input() customClass = '';
  @Input() closeCallback = () => this.show = false;

}
