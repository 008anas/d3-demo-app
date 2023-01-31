import { Component, OnInit, Input } from '@angular/core';

@Component({
  selector: 'sqy-info-popup',
  templateUrl: './info-popup.component.html',
  styleUrls: ['./info-popup.component.scss']
})
export class InfoPopupComponent implements OnInit {

  @Input() content: string = null;

  constructor() { }

  ngOnInit(): void {
  }

}
