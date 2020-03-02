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

class GraphsOptions {
  name: string;
  display: boolean;
  data: any;
  values: string;
  operator: string;
  cutoff?: number;
  color: string;
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
    if (!this._data.results || this._data.results.length) {
      return;
    }
    Object.keys(this._data.results).forEach((key, index) => {
      let value = JSON.parse(this._data.results[key]);
      this.options.push({
        name: key,
        display: true,
        data: value,
        color: this.colors[index],
        values: 'raw',
        operator: this.operators[0].op
      });
    });
  }

  _data: ResultData = null;
  isSearching = false;
  options: GraphsOptions[] = [];
  enzime: string;
  colors = [
    'red',
    'blue',
    'black',
    'green',
    'orange',
    'purple'
  ];
  enzimes = [
    'GAATTC',
    'CTTAAG',
    'CCWGG',
    'GGWCC',
    'GGATCC',
    'CCTAGG',
    'AAGCTT',
    'TTCGAA'
  ];
  operators = [
    { desc: 'Greater than', op: '<' },
    { desc: 'Greater or equal than', op: '<=' },
    { desc: 'Equal', op: '=' },
    { desc: 'Lower than', op: '>' },
    { desc: 'Lower or equal than', op: '>=' }
  ];
  operatorsFn = {
    '<': (a: number, b: number) => a < b,
    '<=': (a: number, b: number) => a <= b,
    '=': (a: number, b: number) => a == b,
    '>': (a: number, b: number) => a > b,
    '>=': (a: number, b: number) => a >= b,
  };

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
      document.querySelectorAll('.matrix-graph').forEach((x: any, i: number) => {
        let data: any;
        data = this.options[i].data.map((v: any) => {
          return {
            pos: v.pos,
            score: this.options[i].values === 'raw' ? v.raw_score : v.norm_score
          }
        });
        x.data = data;
        x.setAttribute('color', this.options[i].color);
      });
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

  setThreshold(value: number, graph_i: number) {
    if (this.options[graph_i] && this.options[graph_i].values && this.options[graph_i].operator && value) {
      let values = this.options[graph_i].data.filter((d: any) => this.operatorsFn[this.options[graph_i].operator](value, this.options[graph_i].values == 'raw' ? d.raw_score : d.norm_score));
      if (values.length) {
        document.getElementById('matrix-graph' + graph_i)['cutoffs'] = values.map((v: { pos: number; }) => v.pos)
      } else {
        this.notify.info('No match was found');
        this.clearThreshold(graph_i);
      }
    }
  }

  changeValues(graph_i: number) {
    const data = this.options[graph_i].data.map((v: any) => {
      return {
        pos: v.pos,
        score: this.options[graph_i].values == 'raw' ? v.raw_score : v.norm_score
      }
    });
    // document.getElementById('matrix-graph' + graph_i)['data'] = data;
  }

  clearThreshold(graph_i: number) {
    document.getElementById('matrix-graph' + graph_i)['cutoffs'] = undefined;
  }

  clearAllThreshold() {
    document.querySelectorAll('.matrix-graph').forEach((x: any) => x.cutoffs = undefined);
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
