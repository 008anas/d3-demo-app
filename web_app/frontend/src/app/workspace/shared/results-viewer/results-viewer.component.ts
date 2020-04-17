import { Component, Input, AfterViewInit, ViewChild, ElementRef } from '@angular/core';
import { ActivatedRoute } from '@angular/router';
import { finalize } from 'rxjs/operators';

import { NzMessageService } from 'ng-zorro-antd/message';
import { NzModalService } from 'ng-zorro-antd/modal';

import { SqrutinyService } from '@services/sqrutiny.service';
import { Track } from 'app/optimizer/shared/track';
import { DisplayValuesComponent } from '../display-values/display-values.component';
import { HistoryService } from '../history.service';
import { FileService } from '@services/file.service';
import { LoaderService } from '@services/loader.service';
import { DisplayThresholdComponent } from '../display-threshold/display-threshold.component';
import { ExportModalComponent } from '../export-modal/export-modal.component';
import { FeatureService } from 'app/optimizer/shared/feature.service';
import { Feature } from 'app/optimizer/shared/feature';
import { SetCutoffComponent } from '../set-cutoff/set-cutoff.component';

class GraphsOptions {
  name?: string;
  alias: string;
  description?: string;
  display: boolean;
  data: any;
  max: number;
  min: number;
  type: string;
  cutoffs?: number[];
  color: string;
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
    this.loader.startLoading();
    this._data.results.forEach((r: any) => {
      this.options.push({
        alias: r.alias,
        display: true,
        data: r.scores.map((d: any) => ({ pos: d.start, score: d.raw_score })),
        max: Math.max.apply(Math, r.scores.map((d: any) => d.raw_score)),
        min: Math.min.apply(Math, r.scores.map((d: any) => d.raw_score)),
        type: 'raw',
        color: this.getRandomColor()
      });
    });
    this.getFeatures();
  }

  _data: any = null; // Original data
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
  @ViewChild('tracksView') tracksElem: ElementRef<HTMLElement>;
  @ViewChild('navigator') nav: ElementRef<HTMLElement>;

  constructor(
    private sqSrvc: SqrutinyService,
    private notify: NzMessageService,
    private modal: NzModalService,
    private fileSrvc: FileService,
    private historySrvc: HistoryService,
    private route: ActivatedRoute,
    private loader: LoaderService,
    private featureSrvc: FeatureService
  ) { }

  ngAfterViewInit() {
    // Init graphs
    if (this._data && this.options) {
      this.nav.nativeElement['length'] = this._data.construct.dna_seq.length;
      document.querySelectorAll('.protvista').forEach((x: any) => x.setAttribute('length', this._data.construct.dna_seq.length));
      this.dnaSeq.nativeElement['data'] = this._data.construct.dna_seq;
      this.proteinSeq.nativeElement['data'] = this.toProteinSeqVw(this._data.construct.protein_seq);
      this.tracksElem.nativeElement['data'] = this.getTrackView(this._data.construct.tracks);
      document.querySelectorAll('.score-graph').forEach((x: any, i: number) => {
        x.data = this.options[i].data;
        x.setAttribute('color', this.options[i].color);
      });
    }
  }

  private getTrackView(tracks: Track[]) {
    return tracks.map((track: Track) => ({
      accession: 'XXXXX',
      start: track.start,
      end: track.end,
      color: track.color
    }));
  }

  private toProteinSeqVw(val: string) {
    return ` ${val.split('').join('  ')} `;
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

  private getRandomColor() {
    const letters = '0123456789ABCDEF';
    let color = '#';
    for (let i = 0; i < 6; i++) {
      color += letters[Math.floor(Math.random() * 16)];
    }
    return color;
  }

  setHighlight(pos: string) {
    if (pos) {
      document.querySelectorAll('.protvista').forEach((x: any) => x.fixedHighlight = pos);
    }
  }

  clearHighlight() {
    document.querySelectorAll('.protvista').forEach((x: any) => x.fixedHighlight = undefined);
  }

  setThreshold(alias: string, value: number, operator: string) {
    if (alias && operator && value !== null) {
      let values = [];
      values = this.getByAlias(this.options, alias).data.filter(d => {
        if (this.operatorsFn[operator](value, d.score)) {
          return {
            pos: d.pos,
            score: d.score
          };
        }
      });
      if (values.length) {
        this.getByAlias(this.options, alias).cutoffs = values;
        document.getElementById(alias)['cutoffs'] = values.map((v: any) => v.pos);
        const element = document.getElementById('cutoffRes' + alias).getElementsByTagName('p')[0];
        element.innerHTML = `<i class="filter icon"></i> ${this.operators.find(o => o.op === operator).desc} ${value}`;
        element.style.visibility = 'visible';
      } else {
        this.clearCutoffs(alias);
        this.notify.info(this.getByAlias(this.options, alias).name + ': No result was found');
      }
    }
  }

  clearCutoffs(graph: string) {
    document.getElementById(graph)['cutoffs'] = undefined;
    this.getByAlias(this.options, graph).cutoffs = null;
    document.getElementById('cutoffRes' + graph).getElementsByTagName('p')[0].innerHTML = '';
    document.getElementById('cutoffRes' + graph).style.visibility = 'hidden';
  }

  clearAllCutoffs() {
    document.querySelectorAll('.score-graph').forEach((x: any) => x.cutoffs = undefined);
    this.options.map(op => op.cutoffs = null);
    document.querySelectorAll('.cutoff').forEach((c: HTMLElement) => c.style.visibility = 'hidden');
  }

  cutoffModal(alias: string) {
    const op = this.getByAlias(this.options, alias);
    this.modal.create({
      nzTitle: 'Set cutoff for ' + op.name,
      nzContent: SetCutoffComponent,
      nzWrapClassName: 'center-modal',
      nzComponentParams: {
        option: op
      },
      nzOnOk: (cmp: SetCutoffComponent) => {
        if (cmp.op !== null && cmp.value !== null) {
          this.setThreshold(alias, cmp.value, cmp.op);
        }
      }
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

  changeColors(alias?: string) {
    if (alias) {
      this.getByAlias(this.options, alias).color = this.getRandomColor();
      document.getElementById(alias).setAttribute('color', this.getByAlias(this.options, alias).color);
    } else {
      this.options.forEach(o => o.color = this.getRandomColor());
      document.querySelectorAll('.score-graph').forEach((x: any, i: number) => x.setAttribute('color', this.options[i].color));
    }
  }

  switchScore(alias: string, type: string) {
    const op = this.getByAlias(this._data.results, alias);
    if (op) {
      op.scores = op.scores.map(s => type === 'raw' ? s.raw_score : s.norm_score);
      document.getElementById(alias)['data'] = op.scores.map((s: any) => ({
        pos: s.start,
        score: type === 'raw' ? s.raw_score : s.norm_score,
        min: Math.min.apply(Math, op.scores),
        max: Math.max.apply(Math, op.scores)
      }));
      this.getByAlias(this.options, alias).type = type;
    }

  }

  getFeatures() {
    this.loader.startLoading();
    this.featureSrvc
      .getAll()
      .subscribe(data => {
        const features = data.map((e: any) => new Feature().deserialize(e));
        this.options.map(o => {
          const ft = features.find(f => f.alias === o.alias);
          if (ft) {
            o.name = ft.name;
            o.description = ft.description;
          }
        });
        this.loader.stopLoading();
      }
      );
  }

  // Export

  exportModal() {
    this.modal.create({
      nzTitle: 'Export wizard',
      nzContent: ExportModalComponent,
      nzComponentParams: {
        items: this.options.map(o => ({ name: o.name, color: o.color, alias: o.alias }))
      },
      nzWidth: 700,
      nzOkText: 'Export',
      nzOnOk: (cmp: ExportModalComponent) => cmp.to === 'gb' ? this.exportToGb(cmp.list.map(l => ({ key: l })) || []) : this.exportToExcel(cmp.list || [])
    });
  }

  exportToGb(options: any) {
    this.loader.startLoading();
    if (options && !Array.isArray(options)) {
      options = [options];
    }
    this.historySrvc.export(this.route.snapshot.data.history.id, { options: options || null })
      .pipe(finalize(() => this.loader.stopLoading()))
      .subscribe((d: any) => {
        this.notify.success('Your download is about to start.');
        this.fileSrvc.saveFileAs(d.data, d.mimetype, d.filename);
      },
        err => this.notify.error(err)
      );
  }

  exportToExcel(alias?: string[]) {
    let data = [];
    if (alias) {
      if (!Array.isArray(alias)) {
        alias = [alias];
      }
      alias.map(a => {
        const op = this.getByAlias(this.options, a);
        data.push({ name: op.name, data: op.data });
      });
    } else {
      data = this._data.results.map(v => ({ name: v.name, data: v.scores }));
    }
    if (data.length > 0) {
      this.fileSrvc.exportAsExcelFile(data, data.length === 1 ? alias[0] : 'scores');
      this.notify.success('Exported! Your download is about to start.');
    } else {
      this.notify.error('Unable to export.');
    }
  }

  private getByAlias(data: any[], alias: string): any {
    return data.find(o => o.alias === alias);
  }
}
