import { Component, Input, AfterViewInit, ViewChild, ElementRef } from '@angular/core';
import { ActivatedRoute } from '@angular/router';
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

class GraphsOptions {
  name: string;
  alias: string;
  display: boolean;
  data: any;
  cutoffs?: number[];
  color: string;
}

class Filter {
  op: string;
  value: number;
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
          }
        }),
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
  fileType = 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet;charset=UTF-8';
  fileExtension = '.xlsx';
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
    { desc: 'Greater than', op: '<' },
    { desc: 'Equal', op: '=' },
    { desc: 'Lower or equal than', op: '>=' },
    { desc: 'Lower than', op: '>' }
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
    private loader: LoaderService
  ) { }

  ngAfterViewInit() {
    this.init();
  }

  init() {
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

  clearHighlight() {
    document.querySelectorAll('.protvista').forEach((x: any) => x.fixedHighlight = undefined);
  }

  setThreshold(graph: string, value: number, operator: string) {
    if (graph && operator) {
      const values = this.options.find(op => op.alias === graph).data.filter((d: any) => this.operatorsFn[operator](value, d.score));
      if (values.length) {
        this.options.find(op => op.alias === graph).cutoffs = values;
        document.getElementById(graph)['cutoffs'] = values.map((v: any) => v.pos);
      } else {
        this.clearThreshold(graph);
        this.notify.info('No match was found');
      }
    }
  }

  clearThreshold(graph: string) {
    document.getElementById(graph)['cutoffs'] = undefined;
    this.options.find(op => op.alias === graph).cutoffs = null;
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
    this.filters.map((f: Filter) => {
      if (f.key && f.op) {
        this.setThreshold(f.key, f.value, f.op);
      }
    });
    this.isVisible = false;
  }

  addFilter() {
    this.filters.push({
      op: this.operators[0].desc,
      key: this.options[0].name,
      value: 0,
      type: 'raw'
    });
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
      this.options.find(o => o.alias === key).color = this.getRandomColor();
      document.getElementById(key).setAttribute('color', this.options.find(o => o.alias === key).color);
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

  export() {
    this.loader.startLoading();
    this.historySrvc.export(this.route.snapshot.data.history.id, { filters: this.filters || null, bulk: this.bulk || false })
      .subscribe((d: any) => {
        this.notify.success('Exported! Your download is about to start.');
        this.fileSrvc.saveFileAs(d.data, d.mimetype, d.filename);
      },
        err => this.notify.error(err),
        () => this.loader.stopLoading()
      );
  }

}
