import { Component, Input } from '@angular/core';

@Component({
  selector: 'sqy-set-cutoff',
  templateUrl: './set-cutoff.component.html',
  styleUrls: ['./set-cutoff.component.scss']
})
export class SetCutoffComponent {

  @Input() option = null;

  operators = [
    { desc: 'Greater or equal than', op: '<=' },
    { desc: 'Greater than', op: '>' },
    { desc: 'Equal', op: '=' },
    { desc: 'Lower or equal than', op: '>=' },
    { desc: 'Lower than', op: '<' }
  ];
  op = this.operators[0].op;
  value: number;

}
