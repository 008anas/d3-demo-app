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
  public colors: string[] = [
    '#000105',
    '#3e6158',
    '#3f7a89',
    '#96c582',
    '#b7d5c4',
    '#bcd6e7',
    '#7c90c1',
    '#9d8594',
    '#dad0d8',
    '#4b4fce',
    '#4e0a77',
    '#a367b5',
    '#ee3e6d',
    '#d63d62',
    '#c6a670',
    '#f46600',
    '#cf0500',
    '#efabbd',
    '#8e0622',
    '#f0b89a',
    '#f0ca68',
    '#62382f',
    '#c97545',
    '#c1800b'
  ];

  constructor(private builder: FormBuilder) {
    this.trackForm = this.builder.group({
      label: [''],
      color: [this.colors[0]],
      sequence: ['', Validators.pattern(Utils.dnaSeqRegex)]
    });
  }

  ngOnChanges(changes: import('@angular/core').SimpleChanges): void {
    if (changes.track && changes.track.currentValue) { this.updateTrackForm(changes.track.currentValue); }
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
    this._track.label = this.label.value;
    this._track.color = this.color.value;
    this._track.sequence = this.sequence.value;
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
