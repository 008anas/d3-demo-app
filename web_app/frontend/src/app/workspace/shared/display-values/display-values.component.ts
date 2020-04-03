import { Component, Input } from '@angular/core';

import { NzMessageService } from 'ng-zorro-antd/message';

import { FileService } from '@services/file.service';

@Component({
  selector: 'sqy-display-values',
  templateUrl: './display-values.component.html',
  styleUrls: ['./display-values.component.scss']
})
export class DisplayValuesComponent {

  _values: number[];
  max: number;
  min: number;

  @Input() set values(values: number[]) {
    this._values = [];
    if (values) {
      this._values = values;
      this.max = Math.max.apply(Math, this._values);
      this.min = Math.min.apply(Math, this._values);
    }
  }

  constructor(
    private fileSrvc: FileService,
    private notify: NzMessageService
  ) { }

  exportToExcel() {
    if (this.values) {
      this.fileSrvc.exportAsExcelFile([{ name: 'Scores', data: this.values }], 'sample');
      this.notify.success('Exported successfully!');
    }
  }

}
