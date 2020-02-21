import { Component, Input, Output, EventEmitter, HostListener, OnDestroy } from '@angular/core';
import { FormGroup, FormBuilder, FormControl } from '@angular/forms';

import { Track } from '../shared/track';
import Utils from 'app/shared/utils';

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
  trackForm = new FormGroup ({
    label: new FormControl(),
    color: new FormControl(),
    sequence: new FormControl()
  });

  constructor() {
    this.dnaRegex = Utils.dnaSeqRegex;
  }

  @HostListener('window:keyup', ['$event'])
  keyEvent(event: KeyboardEvent) {
    if (event.key.toUpperCase() === 'ESCAPE') {
      this.display = false;
    }
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

  onSubmit() {
    console.log(this._track)
    this.onSave.emit(this._track);
    this.toggleSidebar();
  }

  change(pos: number) {
    console.log(this.trackForm)
    this.trackForm.reset();
    // if (form.valid) {
    //   this.onSave.emit(this._track);
    // }
    // if (pos > -1 && pos < this.max + 1) {
    //   this.changePos.emit(pos);
    //   form.reset();
    // }
  }

}
