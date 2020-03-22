import { Component, OnInit, Input, Output, EventEmitter } from '@angular/core';
import { FormGroup, Validators, FormBuilder } from '@angular/forms';

@Component({
  selector: 'sqy-set-highlight',
  templateUrl: './set-highlight.component.html',
  styleUrls: ['./set-highlight.component.scss']
})
export class SetHighlightComponent implements OnInit {

  @Input() max: number;
  @Output() onSet = new EventEmitter<string>();
  highlightForm: FormGroup;

  constructor(private builder: FormBuilder) { }

  get from() {
    return this.highlightForm.get('from');
  }

  get to() {
    return this.highlightForm.get('to');
  }

  ngOnInit(): void {
    this.highlightForm = this.builder.group({
      from: [0, Validators.compose([
        Validators.min(0),
        Validators.required,
        Validators.max(this.max - 1)
      ])],
      to: [this.max, Validators.compose([
        Validators.required,
        Validators.min(1),
        Validators.max(this.max)
      ])]
    });
  }

  set() {
    if (this.from.value < this.to.value) {
      this.onSet.emit(`${this.from.value}:${this.to.value}`);
    }
  }

}
