import { Component, Input } from '@angular/core';

import { NzMessageService } from 'ng-zorro-antd/message';

import { FileService } from '@services/file.service';

@Component({
  selector: 'sqy-display-values',
  templateUrl: './display-values.component.html',
  styleUrls: ['./display-values.component.scss']
})
export class DisplayValuesComponent {

  @Input() headers: string[] = [];
  @Input() values: any[] = [];
  @Input() min: number;
  @Input() max: number;
  @Input() title: string;

  objectKeys = Object.keys;

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
