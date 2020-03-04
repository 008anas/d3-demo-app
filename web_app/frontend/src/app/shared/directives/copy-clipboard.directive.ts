import { Directive, Input, HostListener } from '@angular/core';

import { NzMessageService } from 'ng-zorro-antd/message';

@Directive({
  selector: '[clipboard]'
})
export class CopyClipboardDirective {

  @Input('clipboard')
  public payload: string;

  @Input()
  public context: string;

  constructor(private notify: NzMessageService) { }

  @HostListener('click', ['$event'])
  public onClick(event: MouseEvent): void {
    event.preventDefault();
    if (!this.payload) { return; }

    const listener = (e: ClipboardEvent) => {
      const clipboard = e.clipboardData || window['clipboardData'];
      clipboard.setData('text', this.payload.toString());
      e.preventDefault();
      this.notify.success('Copied!');
    };

    document.addEventListener('copy', listener, false);
    document.execCommand('copy');
    document.removeEventListener('copy', listener, false);
  }

}
