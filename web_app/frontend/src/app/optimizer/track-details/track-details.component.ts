import { Component, Input, Output, EventEmitter, HostListener, OnDestroy } from '@angular/core';
import { Track } from '../shared/track';
import { Base } from '../shared/base';
import Utils from 'src/app/shared/utils';
import { KEY_CODE } from 'src/app/shared/models/key-code';

@Component({
  selector: 'sqy-track-details',
  templateUrl: './track-details.component.html',
  styleUrls: ['./track-details.component.scss']
})
export class TrackDetailsComponent implements OnDestroy {

  @Output() onSave = new EventEmitter<Track>();
  @Output() changePos = new EventEmitter<number>();
  @Input() max: number;
  @Input() set track(x: Track) {
    if (x) {
      this._track = x;
      this.max = this.max || this._track['pos'];
      this.display = true;
    }
  }
  dnaRegex: RegExp;
  display = false;
  _track: Track = null;
  action = 'fix';
  bases: Base[];
  public backgroundColor = '#fff';
  public fontColor = '#222';
  public linkColor = '#4b4fce';

  constructor() {
    this.dnaRegex = Utils.dnaSeqRegex;
  }

  @HostListener('window:keyup', ['$event'])
  keyEvent(event: KeyboardEvent) {
    if (event.keyCode === KEY_CODE.ESCAPE) { this.display = false; }
  }

  toggleSidebar = () => {
    this.display = !this.display;
  }

  ngOnDestroy() {
    this._track = null;
  }

  /*
  * Set color from color picker
  */
  public setColor(color: string) {
    this._track.color = color;
  }

  submit() {
    this.onSave.emit(this._track);
    this.toggleSidebar();
  }

  change(pos: number) {
    if (pos > -1 && pos < this.max + 1) {
      this.changePos.emit(pos);
    }
  }

}
