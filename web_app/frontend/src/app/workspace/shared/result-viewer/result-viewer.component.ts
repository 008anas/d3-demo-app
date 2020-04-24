import { Component, Input, AfterViewInit, ViewChild, ElementRef } from '@angular/core';
import { ActivatedRoute } from '@angular/router';
import { finalize } from 'rxjs/operators';

import { NzMessageService } from 'ng-zorro-antd/message';
import { NzModalService } from 'ng-zorro-antd/modal';
import { ProtvistaSequence } from 'protvista-sequence';
import { ProtvistaSaver } from 'protvista-saver';
import { ProtVistaNavigation } from '@vendor/protvista-navigation';
import { ProtVistaInterproTrack } from 'protvista-interpro-track';

import { SqrutinyService } from '@services/sqrutiny.service';
import { Track } from 'app/optimizer/shared/track';
import { DisplayValuesComponent } from '../display-values/display-values.component';
import { HistoryService } from '../history.service';
import { FileService } from '@services/file.service';
import { LoaderService } from '@services/loader.service';
import { ExportModalComponent } from '../export-modal/export-modal.component';
import { FeatureService } from 'app/optimizer/shared/feature.service';
import { Feature } from 'app/optimizer/shared/feature';
import { SetCutoffComponent } from '../set-cutoff/set-cutoff.component';

class GraphsOption {
  name?: string;
  alias: string;
  description?: string;
  display: boolean;
  data: any;
  max?: number;
  min?: number;
  type: string;
  cutoffs?: number[];
  color: string;
}

class Filter {
  alias: string;
  value: number;
  op: string;
  type: string;
}

@Component({
  selector: 'sqy-result-viewer',
  templateUrl: './result-viewer.component.html',
  styleUrls: ['./result-viewer.component.scss']
})
export class ResultViewerComponent implements AfterViewInit {

  @Input() set data(data: any) {
    this._data = [];
    this.options = [];
    if (!data.result || !data.construct || !data.result.length) {
      return;
    }
    this._data = data;
    this.loader.startLoading();
    this._data.result.map((r: any) => {
      this.options.push({
        alias: r.alias,
        display: true,
        data: r.scores.map((d: any) => ({ pos: d.start, score: d.raw_score })),
        type: 'raw',
        color: this.getRandomColor()
      });
    });
    this.getFeatures();
  }

  _data: any = null; // Original data
  isSearching = false;
  options: GraphsOption[] = [];
  filters: Filter[] = [];
  enzime: string;
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
    { desc: 'Equal to', op: '=' },
    { desc: 'Lower or equal than', op: '>=' },
    { desc: 'Lower than', op: '>' }
  ];
  operatorsFn = {
    '<': (a: number, b: number) => a < b,
    '<=': (a: number, b: number) => a <= b,
    '=': (a: number, b: number) => a === b,
    '>': (a: number, b: number) => a > b,
    '>=': (a: number, b: number) => a >= b,
  };

  @ViewChild('dnaSeq') dnaSeq: ElementRef<ProtvistaSequence>;
  @ViewChild('proteinSeq') proteinSeq: ElementRef<ProtvistaSequence>;
  @ViewChild('tracksView') tracksElem: ElementRef<ProtVistaInterproTrack>;
  @ViewChild('navigator') nav: ElementRef<ProtVistaNavigation>;
  @ViewChild('saver') saver: ElementRef<ProtvistaSaver>;

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
      this.nav.nativeElement.length = this._data.construct.dna_seq.length;
      document.querySelectorAll('.protvista').forEach((x: any) => x.setAttribute('length', this._data.construct.dna_seq.length));
      this.dnaSeq.nativeElement.data = this._data.construct.dna_seq;
      this.proteinSeq.nativeElement.data = this.toProteinSeqVw(this._data.construct.protein_seq);
      this.tracksElem.nativeElement.data = this.getTrackView(this._data.construct.tracks);
      document.querySelectorAll('.score-graph').forEach((x: any, i: number) => {
        x.data = this.options[i].data;
        x.setAttribute('color', this.options[i].color);
      });
      this.saver.nativeElement.preSave = () => {
        this.loader.startLoading();
        const base = document.querySelector('#results-graphs');
        const title = document.createElement('h1');
        title.setAttribute('id', 'tmp_title_element');
        title.innerHTML = 'SQrutiny Result Snapshot';
        base.insertBefore(title, base.firstChild);
        base.querySelectorAll('.no-sv').forEach((element: HTMLElement) => element.style.display = 'none');
        base.querySelectorAll('h3').forEach((element: HTMLElement) => element.style.fontSize = '12px');
      };
      this.saver.nativeElement.postSave = () => {
        const base = document.querySelector('#results-graphs');
        base.removeChild(document.getElementById('tmp_title_element'));
        base.querySelectorAll('.no-sv').forEach((element: HTMLElement) => element.style.display = 'inherit');
        base.querySelectorAll('h3').forEach((element: HTMLElement) => element.style.fontSize = '1.28571429rem');
        this.loader.stopLoading();
      };
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

  private getRandomColor() {
    const letters = '0123456789ABCDEF';
    let color = '#';
    for (let i = 0; i < 6; i++) {
      color += letters[Math.floor(Math.random() * 16)];
    }
    return color;
  }

  displayAll() {
    this.options.map(x => x.display = true);
  }

  hideAll() {
    this.options.map(x => x.display = false);
  }

  setHighlight(pos: string) {
    if (pos) {
      document.querySelectorAll('.protvista').forEach((x: any) => x.fixedHighlight = pos);
    }
  }

  clearHighlight() {
    document.querySelectorAll('.protvista').forEach((x: any) => x.fixedHighlight = undefined);
  }

  clearCutoffs(alias: string) {
    /* tslint:disable:no-string-literal */
    document.getElementById(alias)['cutoffs'] = undefined;
    this.getByAlias(this.options, alias).cutoffs = null;
    this.filters = this.filters.filter(f => f.alias !== alias);
    document.getElementById('cutoff_' + alias).getElementsByTagName('p')[0].innerHTML = '';
    document.getElementById('cutoff_' + alias).style.visibility = 'hidden';
    /* tslint:enable:no-string-literal */
  }

  clearAllCutoffs() {
    document.querySelectorAll('.score-graph').forEach((x: any) => x.cutoffs = undefined);
    this.options.map(op => op.cutoffs = null);
    this.filters = [];
    document.querySelectorAll('.cutoff').forEach((c: HTMLElement) => c.style.visibility = 'hidden');
  }

  cutoffModal(alias: string) {
    const op = this.getByAlias(this.options, alias);
    if (op) {
      this.modal.create({
        nzTitle: 'Set cutoff for ' + op.name,
        nzContent: SetCutoffComponent,
        nzWidth: 450,
        nzWrapClassName: 'center-modal',
        nzComponentParams: {
          option: op
        },
        nzOnOk: (cmp: SetCutoffComponent) => {
          if (cmp.op !== null && cmp.value !== null && !Number.isNaN(cmp.value)) {
            this.setThreshold({ alias, value: cmp.value, op: cmp.op, type: op.type });
          }
        }
      });
    }
  }

  private setThreshold(filter: Filter) {
    if (filter.alias && filter.op && filter.value !== null) {
      let values = [];
      const op = this.getByAlias(this.options, filter.alias);
      values = op.data.filter((d: any) => {
        if (this.operatorsFn[filter.op](filter.value, d.score)) {
          return {
            pos: d.pos,
            score: d.score
          };
        }
      });
      const operator = this.operators.find(o => o.op === filter.op);
      if (values.length) {
        op.cutoffs = values;
        this.filters.push(filter);
        /* tslint:disable:no-string-literal */
        document.getElementById(filter.alias)['cutoffs'] = values.map((v: any) => v.pos);
        /* tslint:enable:no-string-literal */
        const element = document.getElementById('cutoffRes' + filter.alias).getElementsByTagName('p')[0];
        element.innerHTML = `<i class='filter icon'></i> ${operator.desc} ${filter.value}`;
        element.style.visibility = 'visible';
      } else {
        this.clearCutoffs(filter.alias);
        this.notify.info(`No results found for: ${operator.desc} ${filter.value} in ${op.name}`);
      }
    }
  }

  displayScore(alias: string) {
    const data = this.getByAlias(this._data.result, alias);
    this.valuesModal(data.scores, ['Start', 'End', 'Raw values', 'Normalized values'], this.getByAlias(this.options, alias).name + ' data');
  }

  displayCutoffs(alias: string) {
    const op = this.getByAlias(this.options, alias);
    this.valuesModal(op.cutoffs, ['Position', 'Score'], op.name + ' cutoff positions');
  }

  private valuesModal(values: any[], headers: string[], title?: string) {
    this.modal.create({
      nzContent: DisplayValuesComponent,
      nzWrapClassName: 'center-modal',
      nzComponentParams: {
        headers,
        values,
        title
      },
      nzFooter: null
    });
  }

  changeColors(alias?: string) {
    if (alias) {
      const op = this.getByAlias(this.options, alias);
      op.color = this.getRandomColor();
      document.getElementById(alias).setAttribute('color', op.color);
    } else {
      this.options.forEach(o => o.color = this.getRandomColor());
      document.querySelectorAll('.score-graph').forEach((x: any, i: number) => x.setAttribute('color', this.options[i].color));
    }
  }

  switchScore(alias: string, type: string) {
    const data = this.getByAlias(this._data.result, alias);
    if (data) {
      const op = this.getByAlias(this.options, alias);
      this.clearCutoffs(op.alias);
      op.data = data.scores.map((s: any) => ({
        pos: s.start,
        score: type === 'raw' ? s.raw_score : s.norm_score
      }));
      /* tslint:disable:no-string-literal */
      document.getElementById(alias)['data'] = op.data;
      /* tslint:enable:no-string-literal */
      op.type = type;
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

  private getFeatures() {
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
            o.min = ft.genome_min;
            o.max = ft.genome_max;
          }
        });
        this.loader.stopLoading();
      }
      );
  }

  // Export

  exportCutoffToGb(alias: string) {
    const filter = this.getByAlias(this.filters, alias);
    if (filter) {
      this.exportToGb([{ key: alias, filter }]);
    }
  }

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

  exportCutoffToExcel(alias?: string[]) {
    let data = [];
    if (alias) {
      if (!Array.isArray(alias)) {
        alias = [alias];
      }
      alias.map(a => {
        const op = this.getByAlias(this.options, a);
        if (op.cutoffs) {
          data = [{ name: op.name, data: op.cutoffs }];
        }
      });
    } else {
      data = this.options.filter((o: GraphsOption) => !o.cutoffs);
    }
    if (data.length > 0) {
      this.fileSrvc.exportAsExcelFile(data, data.length === 1 ? alias[0] : 'scores');
      this.notify.success('Exported! Your download is about to start.');
    } else {
      this.notify.error('Unable to export.');
    }
  }

  exportToExcel(alias?: string[]) {
    let data = [];
    if (alias) {
      if (!Array.isArray(alias)) {
        alias = [alias];
      }
      alias.map(a => {
        const op = this.getByAlias(this._data.result, a);
        data.push({ name: op.alias, data: op.scores });
      });
    } else {
      data = this._data.result.map((v: any) => ({ name: v.alias, data: v.scores }));
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
