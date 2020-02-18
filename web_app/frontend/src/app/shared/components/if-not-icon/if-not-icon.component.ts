import { Component, Input } from '@angular/core';

@Component({
  selector: 'sqy-if-not-icon',
  template: '{{value}} <i class="{{icon}} {{iconColor}} icon" [hidden]="value && value != 0"></i>'
})
export class IfNotIconComponent {

  @Input() value: any;
  @Input() icon = 'ban';
  @Input() iconColor = 'red';
}
