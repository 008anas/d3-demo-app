import { Component, Input } from '@angular/core';

@Component({
  selector: 'sqy-export-modal',
  templateUrl: './export-modal.component.html',
  styleUrls: ['./export-modal.component.scss']
})
export class ExportModalComponent {

  list: any[] = [];
  _items: any[];
  to = 'gb';

  @Input() set items(items: any[]) {
    if (items.length > 0) {
      this._items = items;
      this.list = this._items.map(f => f.alias);
    }
  }

  featChange(alias: string, isChecked: boolean) {
    isChecked ? this.list.push(alias) : this.list = this.list.filter(f => f !== alias);
  }

}
