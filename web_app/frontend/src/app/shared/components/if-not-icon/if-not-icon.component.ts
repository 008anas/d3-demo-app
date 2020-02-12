import { Component, Input } from '@angular/core';

@Component({
  selector: 'sqy-if-not-icon',
  templateUrl: './if-not-icon.component.html',
  styleUrls: ['./if-not-icon.component.scss']
})
export class IfNotIconComponent {

  @Input() value: any;
  @Input() icon = 'ban';
  @Input() iconColor = 'red';

}
