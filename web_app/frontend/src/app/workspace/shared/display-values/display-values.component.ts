import { Component, Input } from '@angular/core';

import { NzMessageService } from 'ng-zorro-antd/message';

import { FileService } from '@services/file.service';

@Component({
  selector: 'sqy-display-values',
  templateUrl: './display-values.component.html',
  styleUrls: ['./display-values.component.scss']
})
export class DisplayValuesComponent {

  _values: any[];
  max_raw: number;
  min_raw: number;
  max_norm: number;
  min_norm: number;

  @Input() set values(values: any[]) {
    this._values = [];
    if (values) {
      this._values = values;
      const raw = this._values.map(v => v.raw_score);
      const norm = this._values.map(v => v.norm_score);
      this.max_raw = Math.max.apply(Math, raw);
      this.min_raw = Math.min.apply(Math, raw);
      this.max_norm = Math.max.apply(Math, norm);
      this.min_norm = Math.min.apply(Math, norm);
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
