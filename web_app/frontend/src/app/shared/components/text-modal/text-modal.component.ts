import { Component, Input } from '@angular/core';

import { NotifyService } from '@services/notify.service';

@Component({
  selector: 'sqy-text-modal',
  templateUrl: './text-modal.component.html',
  styleUrls: ['./text-modal.component.scss']
})
export class TextModalComponent {

  @Input() txt: string = null;
  @Input() description: string = null;
  @Input() copy = true;

  constructor(
    private notify: NotifyService
  ) { }

  notifySuccess(msg: string) {
    this.notify.success(msg, 'top-right');
  }

}
