import { Component, Input, AfterViewInit, ViewChild, ElementRef } from '@angular/core';
import { ActivatedRoute } from '@angular/router';
import { FormGroup, Validators, FormBuilder } from '@angular/forms';
import { finalize } from 'rxjs/operators';

import { NzMessageService } from 'ng-zorro-antd/message';
import { NzModalService } from 'ng-zorro-antd/modal';

import { SqrutinyService } from '@services/sqrutiny.service';
import Utils from 'app/shared/utils';
import { Track } from 'app/optimizer/shared/track';
import { DisplayValuesComponent } from '../display-values/display-values.component';
import { HistoryService } from '../history.service';
import { FileService } from '@services/file.service';
import { LoaderService } from '@services/loader.service';
import { DisplayThresholdComponent } from '../display-threshold/display-threshold.component';

class GraphsOptions {
  name: string;
  alias: string;
  display: boolean;
  data: any;
  type: string;
  cutoffs?: number[];
  color: string;
}

class Filter {
  op: string;
  value?: number;
  key: string;
  type: string;
}

@Component({
  selector: 'sqy-results-viewer',
  templateUrl: './results-viewer.component.html',
  styleUrls: ['./results-viewer.component.scss']
})
export class ResultsViewerComponent implements AfterViewInit {

  @Input() set data(data: any) {
    this._data = data;
    this.options = [];
    if (!this._data.results || !this._data.construct || !this._data.results.length) {
      return;
    }
    this._data.results.forEach((r) => {
      this.options.push({
        name: r.name,
        alias: r.alias,
        display: true,
        data: r.scores.map((d: any) => {
          return {
            pos: d.start,
            score: d.raw_score
          };
        }),
        type: 'raw',
        color: this.getRandomColor()
      });
    });
  }

  isVisible = false;
  _data: any = null; // Original data
  isSearching = false;
  options: GraphsOptions[] = [];
  enzime: string;
  filters: Filter[] = [];
  bulk = false;
  highlightForm: FormGroup;
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
    { desc: 'Greater or equal than', op: '<=' },
    { desc: 'Greater than', op: '>' },
    { desc: 'Equal', op: '=' },
    { desc: 'Lower or equal than', op: '>=' },
    { desc: 'Lower than', op: '<' }
  ];
  operatorsFn = {
    '<': (a: number, b: number) => a < b,
    '<=': (a: number, b: number) => a <= b,
    '=': (a: number, b: number) => a == b,
    '>': (a: number, b: number) => a > b,
    '>=': (a: number, b: number) => a >= b,
  };

  @ViewChild('dnaSeq') dnaSeq: ElementRef<HTMLElement>;
  @ViewChild('proteinSeq') proteinSeq: ElementRef<HTMLElement>;
  @ViewChild('tracksView') tracksComponent: ElementRef<HTMLElement>;

  constructor(
    private sqSrvc: SqrutinyService,
    private notify: NzMessageService,
    private modal: NzModalService,
    private fileSrvc: FileService,
    private historySrvc: HistoryService,
    private route: ActivatedRoute,
    private loader: LoaderService,
    private builder: FormBuilder
  ) { }

  ngAfterViewInit() {
    this.highlightForm = this.builder.group({
      from: [0, Validators.required],
      to: [this._data.construct.dna_seq.length, Validators.required]
    });
    this.init();
  }

  private init() {
    if (this._data && this.options) {
      document.querySelectorAll('.protvista').forEach((x: any) => x.setAttribute('length', this._data.construct.dna_seq.length));
      this.dnaSeq.nativeElement['data'] = this._data.construct.dna_seq;
      this.proteinSeq.nativeElement['data'] = ResultsViewerComponent.proteinSeqProtvista(this._data.construct.protein_seq);
      this.tracksComponent.nativeElement['data'] = ResultsViewerComponent.getTrackView(this._data.construct.tracks);
      document.querySelectorAll('.score-graph').forEach((x: any, i: number) => {
        x.data = this.options[i].data;
        x.setAttribute('color', this.options[i].color);
      });
      this.dnaSeq.nativeElement.click(); // Fix init issue. Simulate click to show data.
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

  setHighlight() {
    document.querySelectorAll('.protvista').forEach((x: any) => x.fixedHighlight = `${this.highlightForm.value.from}:${this.highlightForm.value.to}`);
  }

  clearHighlight() {
    document.querySelectorAll('.protvista').forEach((x: any) => x.fixedHighlight = undefined);
  }

  setThreshold(graph: string, value: number, operator: string, type: string) {
    if (graph && operator) {
      const values = [];
      this._data.results.find((r: any) => r.alias === graph).scores.map((s: any) => {
        const score = type === 'raw' ? s.raw_score : s.norm_score;
        if (this.operatorsFn[operator](value, score)) {
          values.push({ pos: s.start, score });
        }
      });
      if (values.length) {
        this.options.find(op => op.alias === graph).cutoffs = values;
        document.getElementById(graph)['cutoffs'] = values.map((v: any) => v.pos);
        const element = document.getElementById('cutoffRes' + graph).getElementsByTagName('p')[0]
        element.innerHTML = `<i class="filter icon"></i> ${this.operators.find(o => o.op === operator).desc} ${value}`;
        element.style.visibility = 'visible';
      } else {
        this.clearThreshold(graph);
        this.notify.info(this.getByAlias(graph).name + ': No result was found');
      }
    }
  }

  clearThreshold(graph: string) {
    document.getElementById(graph)['cutoffs'] = undefined;
    this.getByAlias(graph).cutoffs = null;
    this.filters = this.filters.filter(f => f.key !== graph);
    document.getElementById('cutoffRes' + graph).getElementsByTagName('p')[0].innerHTML = '';
    document.getElementById('cutoffRes' + graph).style.visibility = 'hidden';
  }

  clearAllThreshold() {
    document.querySelectorAll('.score-graph').forEach((x: any) => x.cutoffs = undefined);
    this.options.map(op => op.cutoffs = null);
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

  applyFilters() {
    this.clearAllThreshold();
    this.filters.map((f: Filter) => {
      if (f.key && f.op) {
        this.setThreshold(f.key, f.value, f.op, f.type);
      }
    });
    this.isVisible = false;
  }

  addFilter() {
    if (this.filters.length < this.options.length) {
      this.filters.push({
        op: this.operators[0].op,
        key: this.options[this.filters.length].alias,
        type: 'raw'
      });
    }
  }

  valuesModal(key: string) {
    this.modal.create({
      nzContent: DisplayValuesComponent,
      nzWrapClassName: 'center-modal',
      nzComponentParams: {
        values: this._data.results.find(r => r.alias === key).scores
      },
      nzFooter: null
    });
  }

  thresholdModal(key: string) {
    this.modal.create({
      nzContent: DisplayThresholdComponent,
      nzWrapClassName: 'center-modal',
      nzComponentParams: {
        values: this.options.find(r => r.alias === key).cutoffs
      },
      nzFooter: null
    });
  }

  exportToExcel(key?: string) {
    if (key) {
      this.fileSrvc.exportAsExcelFile(this._data.results.find(r => r.alias === key).scores, 'sample');
    } else {
      // TODO: export all
    }
    this.notify.success('Exported successfully!');
  }

  changeColors(key?: string) {
    if (key) {
      this.getByAlias(key).color = this.getRandomColor();
      document.getElementById(key).setAttribute('color', this.getByAlias(key).color);
    } else {
      this.options.forEach(o => o.color = this.getRandomColor());
      document.querySelectorAll('.score-graph').forEach((x: any, i: number) => x.setAttribute('color', this.options[i].color));
    }
  }

  getRandomColor() {
    const letters = '0123456789ABCDEF';
    var color = '#';
    for (var i = 0; i < 6; i++) {
      color += letters[Math.floor(Math.random() * 16)];
    }
    return color;
  }

  getByAlias(alias: string): GraphsOptions {
    return this.options.find(o => o.alias === alias);
  }

  exportThreshold(key?: string) {
    if (key) {
      const f = this.filters.filter(f => f.key === key);
      this.export(f);
    }
  }

  export(filters: Filter[]) {
    this.loader.startLoading();
    this.historySrvc.export(this.route.snapshot.data.history.id, { filters: filters || null, bulk: this.bulk || false })
      .pipe(finalize(() => this.loader.stopLoading()))
      .subscribe((d: any) => {
        this.notify.success('Your download is about to start.');
        this.fileSrvc.saveFileAs(d.data, d.mimetype, d.filename);
      },
        err => this.notify.error(err)
      );
  }

  switchValues(type: string, alias: string) {
    document.getElementById(alias)['data'] = this._data.results.find(r => r.alias === alias).scores.map((s: any) => {
      return {
        pos: s.start,
        score: type === 'raw' ? s.raw_score : s.norm_score
      }
    });
    this.options.find(o => o.alias === alias).type = type;
  }

}
