import { Component, Input } from '@angular/core';

import { Track } from 'app/optimizer/shared/track';

@Component({
  selector: 'sqy-track-display',
  templateUrl: './track-display.component.html',
  styleUrls: ['./track-display.component.scss']
})
export class TrackDisplayComponent {

  @Input() track: Track;
}
