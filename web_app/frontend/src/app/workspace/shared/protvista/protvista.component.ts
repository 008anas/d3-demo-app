import { Component, Input, AfterViewInit, ViewChild, ElementRef } from '@angular/core';
import { finalize } from 'rxjs/operators';

import { NzMessageService } from 'ng-zorro-antd/message';

import { Construct } from '@models/construct';
import { Track } from 'app/optimizer/shared/track';
import { SqrutinyService } from '@services/sqrutiny.service';
import Utils from 'app/shared/utils';

class ResultData {
  construct: Construct;
  results: any;
}

@Component({
  selector: 'sqy-protvista',
  templateUrl: './protvista.component.html',
  styleUrls: ['./protvista.component.scss']
})
export class ProtvistaComponent implements AfterViewInit {

  @Input() set data(data: ResultData) {
    this._data = data;
    this.options = [];
    for (const key in this._data.results) {
      const value = this._data.results[key];
      if (value) {
        this.options.push({
          name: key,
          display: true,
          data: value
        });
      }
    }
  }

  _data: ResultData = null;
  isSearching = false;
  options: any[] = [];
  enzimes = ['GAATTC', 'CTTAAG', 'CCWGG', 'GGWCC', 'GGATCC', 'CCTAGG', 'AAGCTT', 'TTCGAA'];
  enzime: string;

  @ViewChild('dnaSeq') manager: ElementRef<HTMLElement>;

  constructor(
    private sqSrvc: SqrutinyService,
    private notify: NzMessageService
  ) { }

  ngAfterViewInit() {
    this.init();
  }

  init() {
    if (this._data) {
      document.querySelectorAll('.protvista').forEach((x: any) => x.setAttribute('length', this._data.construct.dna_seq.length));
      document.querySelector('#dna_sequence')['data'] = this._data.construct.dna_seq;
      document.querySelector('#protein_sequence')['data'] = ProtvistaComponent.proteinSeqProtvista(this._data.construct.protein_seq);
      document.querySelector('#interpro-track-residues')['data'] = ProtvistaComponent.getTrackView(this._data.construct.tracks);
      document.querySelectorAll('.var-graph').forEach((x: any, i: number) => x.data = JSON.parse(this.options[i].data));
      this.manager.nativeElement.click();
    }
  }

  searchMotif(motif: any) {
    if (motif) {
      this.isSearching = true;
      this.sqSrvc.motifInSeq(this._data.construct.dna_seq, motif)
        .pipe(finalize(() => this.isSearching = false))
        .subscribe(
          (data: any) => {
            if (data.count > 0) {
              let positions = '';
              data.data.map((pos: number[]) => positions += `,${pos[0]}:${pos[1]}`);
              if (positions.charAt(0) === ',') {
                positions = positions.substring(1);
              }
              document.querySelectorAll('.protvista').forEach((x: any) => x.fixedHighlight = positions);
            } else {
              this.notify.error('No match was found');
            }
          },
          err => this.notify.error(err)
        );
    }
  }

  displayAll() {
    this.options.map(x => x.display = true);
  }

  hideAll() {
    this.options.map(x => x.display = false);
  }

  move(x: number, i: number) {
    const pos = x + i;
    if (-1 < pos && pos <= this.options.length - 1) {
      this.options = Utils.array_move(Object.assign([], this.options), x, pos);
    }
  }

  clearHighlight() {
    document.querySelectorAll('.protvista').forEach((x: any) => x.fixedHighlight = undefined);
  }

  static getTrackView(tracks: Track[]) {
    return tracks.map((track: Track) => {
      {
        return {
          accession: 'XXXXX',
          start: track.start,
          end: track.end,
          color: track.color
        };
      }
    });
  }

  static proteinSeqProtvista(val: string) {
    return ` ${val.split('').join('  ')} `;
  }

}
