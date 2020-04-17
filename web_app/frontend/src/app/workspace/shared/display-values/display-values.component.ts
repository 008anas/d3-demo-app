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
  maxRaw: number;
  minRaw: number;
  maxNorm: number;
  minNorm: number;

  @Input() set values(values: any[]) {
    this._values = [];
    if (values) {
      this._values = values;
      const raw = this._values.map(v => v.raw_score);
      const norm = this._values.map(v => v.norm_score);
      this.maxRaw = Math.max.apply(Math, raw);
      this.minRaw = Math.min.apply(Math, raw);
      this.maxNorm = Math.max.apply(Math, norm);
      this.minNorm = Math.min.apply(Math, norm);
    }
  }

  constructor(
    private fileSrvc: FileService,
    private notify: NzMessageService
  ) { }

  exportToExcel() {
    if (this.values) {
      this.fileSrvc.exportAsExcelFile([{ name: 'Scores', data: this.values }]);
      this.notify.success('Exported successfully!');
    }
  }

}
