import { Component, Input, Output, EventEmitter, OnDestroy, OnChanges } from '@angular/core';
import { FormGroup, FormBuilder, Validators } from '@angular/forms';

import { Track } from '../shared/track';
import Utils from 'app/shared/utils';

@Component({
  selector: 'sqy-track-details',
  templateUrl: './track-details.component.html',
  styleUrls: ['./track-details.component.scss']
})
export class TrackDetailsComponent implements OnChanges, OnDestroy {

  @Output() onSave = new EventEmitter<Track>();
  @Output() changePos = new EventEmitter<number>();
  @Input() max: number;
  @Input() set track(x: Track) {
    if (x) {
      this._track = x;
      this.max = this.max || this._track['pos'];
      this.updateTrackForm(x);
      this.display = true;
    }
  }
  display = false;
  _track: Track = null;
  trackForm: FormGroup;

  constructor(private builder: FormBuilder) {
    this.trackForm = this.builder.group({
      label: [''],
      color: [''],
      sequence: ['', Validators.pattern(Utils.dnaSeqRegex)]
    });
  }

  ngOnChanges(changes: import('@angular/core').SimpleChanges): void {
    if (changes.track.currentValue) { this.updateTrackForm(changes.track.currentValue); }
  }

  ngOnDestroy() {
    this._track = null;
  }

  get label() {
    return this.trackForm.get('label');
  }

  get color() {
    return this.trackForm.get('color');
  }

  get sequence() {
    return this.trackForm.get('sequence');
  }

  updateTrackForm(track: Track) {
    this.trackForm.patchValue({
      label: track.label,
      color: track.color,
      sequence: track.sequence
    });
  }

  /*
  * Set color from color picker
  */
  setColor(color: string) {
    this.trackForm.value.color = color;
  }

  onSubmit() {
    this.onSave.emit(this._track);
    this.display = false;
  }

  updateTrack() {
    this._track.label = this.trackForm.value.label;
    this._track.color = this.trackForm.value.color;
    this._track.sequence = this.trackForm.value.sequence;
    this.onSave.emit(this._track);
  }

  change(pos: number) {
    if (this.trackForm.touched && this.trackForm.valid) {
      this.updateTrack();
    }
    if (pos > -1 && pos < this.max + 1) {
      this.changePos.emit(pos);
      this.trackForm.reset();
    }
  }

}
