import { Component, OnInit, Input } from '@angular/core';

const MAX = 100;
const MIN = 0;

@Component({
  selector: 'sqy-progress-bar',
  templateUrl: './progress-bar.component.html',
  styleUrls: ['./progress-bar.component.scss']
})
export class ProgressBarComponent implements OnInit {

  @Input() progress = 0;
  @Input() msg = 'Funded';
  max = MAX;
  min = MIN;

  constructor() { }

  ngOnInit() {
  }

  increment() {
    if (this.progress < MAX) this.progress++;
  }

  decrement() {
    if (this.progress > MIN) this.progress--;
  }

}
