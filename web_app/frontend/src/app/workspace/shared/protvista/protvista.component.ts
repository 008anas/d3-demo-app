import { Component, OnInit, Input, Output, EventEmitter } from '@angular/core';
import { Construct } from 'src/app/construct/shared/construct';

@Component({
  selector: 'sqy-protvista',
  templateUrl: './protvista.component.html',
  styleUrls: ['./protvista.component.scss']
})
export class ProtvistaComponent implements OnInit {

  @Input() set data(data: Construct) {
    if (data) {
      this.construct = data;
      this.isLoading = true;
      this.initGraph();
    }
  }
  @Input() loading = false;
  @Input() status = 'active';
  @Output() onRefresh = new EventEmitter<void>();
  isLoading = false;
  isSearching = false;
  construct: Construct = null;

  constructor() { }

  ngOnInit() {
  }

  searchMotif(event: any) {
    if (event.target.value) {
      this.isSearching = true;
    }
  }

  clearHighlight() {
    document.querySelectorAll('.protvista').forEach((x: any) => x.fixedHighlight = undefined);
  }

  refresh() {
    this.onRefresh.emit();
  }

  initGraph() {
  }

}
