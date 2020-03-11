import { Component, Input } from '@angular/core';

import { NzMessageService } from 'ng-zorro-antd/message';

import { FileService } from '@services/file.service';

@Component({
  selector: 'sqy-display-values',
  templateUrl: './display-values.component.html',
  styleUrls: ['./display-values.component.scss']
})
export class DisplayValuesComponent {

  constructor(
    private fileSrvc: FileService,
    private notify: NzMessageService
  ) { }

  @Input() values: any[];

  exportToExcel() {
    if (this.values) {
      this.fileSrvc.exportAsExcelFile(this.values, 'sample');
      this.notify.success('Exported successfully!');
    }
  }

}
