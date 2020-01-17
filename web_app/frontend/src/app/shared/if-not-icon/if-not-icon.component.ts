import { Component, Input } from '@angular/core';

@Component({
  selector: 'sqy-if-not-icon',
  templateUrl: './if-not-icon.component.html',
  styleUrls: ['./if-not-icon.component.scss']
})
export class IfNotIconComponent {

  @Input() value: any;
  @Input() icon: string = 'ban';
  @Input() iconColor: string = 'red';

}
