import { Component, Input } from '@angular/core';

import { NzMessageService } from 'ng-zorro-antd/message';

import { ExcelService } from '@services/excel.service';

@Component({
  selector: 'sqy-display-values',
  templateUrl: './display-values.component.html',
  styleUrls: ['./display-values.component.scss']
})
export class DisplayValuesComponent {

  constructor(
    private excelSrvc: ExcelService,
    private notify: NzMessageService
  ) { }

  @Input() values: any[];

  exportToExcel() {
    if (this.values) {
      this.excelSrvc.exportAsExcelFile(this.values, 'sample');
      this.notify.success('Exported successfully!');
    }
  }

}
