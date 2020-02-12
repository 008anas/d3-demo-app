import { Component, Input } from '@angular/core';

import { Errors } from '../../models/error';

@Component({
  selector: 'sqy-errors-list',
  templateUrl: './errors-list.component.html',
  styleUrls: ['./errors-list.component.scss']
})
export class ErrorsListComponent {

  formattedErrors: Array<string> = [];

  @Input()
  set errors(errorList: Errors) {
    this.formattedErrors = Object.keys(errorList || {})
      .map(key => errorList[key]);
  }

  get errorList() { return this.formattedErrors; }

}
