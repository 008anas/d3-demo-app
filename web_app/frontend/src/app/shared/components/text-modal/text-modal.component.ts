import { Component, Input } from '@angular/core';

@Component({
  selector: 'sqy-text-modal',
  templateUrl: './text-modal.component.html',
  styleUrls: ['./text-modal.component.scss']
})
export class TextModalComponent {

  @Input() txt: string = null;
  @Input() description: string = null;
  @Input() copy = true;

}
