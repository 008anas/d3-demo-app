import { Component, OnInit, HostListener, OnDestroy } from '@angular/core';
import { Router, ActivatedRoute } from '@angular/router';
import { Subscription } from 'rxjs/internal/Subscription';
import { HttpClient } from '@angular/common/http';
import { finalize } from 'rxjs/operators';

import { saveAs } from 'file-saver';

import { environment as env } from 'src/environments/environment';
import { SpecieService } from 'src/app/shared/services/specie.service';
import { KEY_CODE } from 'src/app/shared/models/key-code';
import { Specie } from 'src/app/shared/models/specie';
import { Construct } from 'src/app/construct/shared/construct';
import { Track } from '../shared/track';
import { TrackService } from '../shared/track.service';
import { NotifyService } from 'src/app/shared/services/notify.service';
import { ConstructService } from 'src/app/construct/shared/construct.service';
import { UserHistory } from 'src/app/workspace/shared/user-history';
import Utils from 'src/app/shared/utils';

class Category {
  name: string;
  elements: Track[];
}

@Component({
  selector: 'sqy-sketcher',
  templateUrl: './sketcher.component.html',
  styleUrls: ['./sketcher.component.scss']
})
export class SketcherComponent implements OnInit, OnDestroy {

  sub: Subscription;
  tracks: Track[] = [];
  track: Track = null;
  history: UserHistory = new UserHistory();
  trackHovered: Track = null;
  hoveredName: string = null;
  isLoading = false;
  response: any = null;
  categories: Category[] = [];
  specie: Specie = new Specie();
  species: Specie[] = [];
  isSubmitting = false;
  submitted = false;
  showPicker = false;
  zoom = 75;
  isTracksLoading = false;
  construct: Construct = new Construct();
  showIndexes = true;
  view = 'general';
  locked = false;
  search: string;

  @HostListener('window:keyup', ['$event'])
  keyEvent(event: KeyboardEvent) {
    if (event.keyCode === KEY_CODE.ESCAPE) {
      this.showPicker = false;
    }
  }

  constructor(
    private route: ActivatedRoute,
    private specieSrvc: SpecieService,
    private trackSrvc: TrackService,
    private constructSrvc: ConstructService,
    private http: HttpClient,
    private router: Router,
    private notify: NotifyService
  ) {
    this.construct.tracks = [];
    this.construct.sequence = '';
  }

  ngOnInit() {
    this.isLoading = true;
    this.sub = this.route.queryParams.subscribe(params => {
      this.construct.id = params.id || null;
      this.specie.slug = params.specie || null;
    });
    if (this.specie.slug) {
      this.getSpecie();
    }
    if (this.construct.id) {
      this.getConstruct();
    } else {
      this.new();
    }
    this.getTracks();
    this.getSpecies();
  }

  ngOnDestroy() {
    if (this.sub) { this.sub.unsubscribe(); }
  }

  getSpecie() {
    this.isLoading = true;
    this.specieSrvc.getBySlug(this.specie)
      .pipe(finalize(() => this.isLoading = false))
      .subscribe(data => this.specie.deserialize(data));
  }

  getConstruct() {
    if (this.construct.id) {
      this.isLoading = true;
      this.constructSrvc.getById(this.construct.id)
        .pipe(finalize(() => this.isLoading = false))
        .subscribe(data => this.construct = new Construct().deserialize(data));
    }
  }

  getExampleConstruct() {
    this.isLoading = true;
    this.constructSrvc.getExample()
      .pipe(finalize(() => this.isLoading = false))
      .subscribe(data => data.length ? this.construct = new Construct().deserialize(data[0]) : this.notify.warn('Sorry but no example construct was found'),
        err => this.notify.warn(err.msg || 'Unable to load model construct'));
  }

  getSpecies() {
    this.isLoading = true;
    this.specieSrvc.getAll()
      .pipe(finalize(() => this.isLoading = false))
      .subscribe(
        data =>
          this.species = data.map((e: any) => {
            const specie = new Specie().deserialize(e);
            if (specie.default && !this.specie.slug) {
              this.specie = Object.assign({}, specie);
            }
            return specie;
          })
      );
  }

  getTracksByCategories() {
    this.isLoading = true;
    this.trackSrvc.getByCategories()
      .pipe(finalize(() => this.isLoading = false))
      .subscribe(data => this.categories = data);
  }

  getTracks() {
    this.isLoading = true;
    this.trackSrvc.getAll()
      .pipe(finalize(() => this.isLoading = false))
      .subscribe(data => {
        this.tracks = data;
        this.new();
      });
  }

  new() {
    if (this.tracks.some(e => e.default)) {
      this.construct.tracks = [Object.assign({}, this.tracks.find(e => e.default))];
    }
  }

  clear() {
    if (confirm('Are you sure you want to clear all tracks from construct? These action cannot be reverted.')) {
      this.construct.tracks = [];
      this.construct.sequence = '';
    }
  }

  removeTrack(x: Track) {
    const i: number = this.construct.tracks.indexOf(x);
    if (i !== -1) {
      this.construct.tracks.splice(i, 1);
    }
  }

  someSelected() {
    return this.categories.some(c => c.elements.some(t => t.selected));
  }

  addTracks() {
    if (this.someSelected() && !this.locked) {
      this.categories.map(c => {
        c.elements.map(e => {
          if (e.selected) {
            e.selected = false;
            this.construct.tracks.push(Object.assign({}, e));
          }
        });
      });
      this.showPicker = false;
    }
  }

  exampleConstruct() {
    let flag = true;
    if (this.construct.tracks.length > 0) {
      if (!confirm('You really want to load an example construct? You\'re gonna lose all actual data. Proceed?')) {
        flag = false;
        // this.construct.name = 'Example construct';
        // this.construct.tracks = Object.assign([], this.tracks.slice(0, 4));
        // this.construct.tracks[0].label = 'First track';
        // this.construct.tracks[0].sequence = 'CGGCGAGCGGCGAGTAACGGCGAGCGGCGAGTAAAATATATAAAATGAGCGGAGAGCGCGCGGCGAGCGGCGAGTAACGGCGAGCGGCGAGTAAAATATATAAAATGAGCGGAGAGCGCG';
        // this.construct.tracks[0].color = '#96c582';
        // this.construct.sequence += this.construct.tracks[0].sequence;
        // this.construct.tracks[1].label = 'Second track';
        // this.construct.tracks[1].sequence = 'CGGCGAGCGGCGAGTAACGGCGAGCGGCGAGTAAAATATATAAAATGAGCGGAGAGCGCGCGGCGAGCGGCGAGTAACGGCGAGCGGCGAGTAAAATATATAAAATGAGCGGAGAGCGCG';
        // this.construct.tracks[1].color = '#cf0500';
        // this.construct.sequence += this.construct.tracks[1].sequence;
        // this.construct.tracks[2].label = 'Third track';
        // this.construct.tracks[2].sequence = 'CGGCGAGCGGCGAGTAACGGCGAGCGGCGAGTAAAATATATAAAATGAGCGGAGAGCGCGCGGCGAGCGGCGAGTAACGGCGAGCGGCGAGTAAAATATATAAAATGAGCGGAGAGCGCG';
        // this.construct.tracks[2].color = '#4b4fce';
        // this.construct.sequence += this.construct.tracks[2].sequence;
        // this.construct.tracks[3].label = 'Fourth track';
        // this.construct.tracks[3].sequence = 'CGGCGAGCGGCGAGTAACGGCGAGCGGCGAGTAAAATATATAAAATGAGCGGAGAGCGCGCGGCGAGCGGCGAGTAACGGCGAGCGGCGAGTAAAATATATAAAATGAGCGGAGAGCGCG';
        // this.construct.tracks[3].color = '#f0ca68';
        // this.construct.sequence += this.construct.tracks[3].sequence;
      }
    }

    if (flag) { this.getExampleConstruct(); }
  }

  submit() {
    if (this.checkTracks()) {
      this.response = null;
      this.isSubmitting = true;
      this.construct['specie_tax_id'] = this.specie.tax_id;
      this.http.post(`${env.endpoints.api}/optimize_seq/from-sketch`, this.construct)
        .pipe(finalize(() => this.isSubmitting = false))
        .subscribe(
          (data: UserHistory) => {
            this.history = new UserHistory().deserialize(data);
            this.submitted = true;
            setTimeout(() => {
              this.submitted = false;
              this.router.navigate(['/workspace', this.history.id]);
            }, 3000);
          },
          err => {
            this.notify.error(err.msg, 'bottom-right');
          });
    }
  }

  // Track Details Sidebar

  openSidebar(e: Track, i: number) {
    if (!this.locked) {
      this.track = Object.assign({}, e);
      this.track['pos'] = i;
    }
  }

  addTrack(track: Track) {
    if (track['pos'] > -1) {
      this.construct.tracks[track['pos']] = track;
      this.construct.tracks[track['pos']].start = this.construct.sequence.length + 1;
      this.construct.tracks[track['pos']].end = this.construct.sequence.length + track.sequence.length;
      this.construct.sequence += track.sequence;
      delete this.construct.tracks[track['pos']]['invalid']; // Now is valid
    }
    this.track = null;
  }

  changeTrack(pos: number) {
    if (pos > -1 && this.construct.tracks[pos]) {
      this.openSidebar(this.construct.tracks[pos], pos);
    }
  }

  moveTrack(x: number, i: number) {
    const pos = x + i;
    if (-1 < pos && pos <= this.construct.tracks.length - 1) { this.construct.tracks = Utils.array_move(this.construct.tracks, x, pos); }
  }

  checkTracks() {
    let flag = true;
    this.construct.tracks.map(t => {
      if (!t.sequence) {
        t['invalid'] = true;
        flag = false;
      }
    });
    return flag;
  }

  toggleSelection(event: { target: { checked: boolean; }; }) {
    this.categories.map(c => c.elements.map(e => e.selected = event.target.checked));
  }

  // Export / Save Construct

  downloadAs(op: string) {
    if (this.construct.tracks.length) {

      let data: BlobPart;
      let ext: string;

      switch (op.toUpperCase()) {
        case 'GENBANK':
          data = Utils.jsonToGenbank(this.construct);
          ext = 'gbk';
          break;
        case 'FASTA':
          data = Utils.jsonToFasta([this.construct]);
          ext = 'fasta';
          break;
        case 'XLSX':
          // TODO:
          break;
        case 'JSON':
          data = JSON.stringify({ construct: this.construct });
          ext = 'json';
          break;
      }

      if (data) {
        const blob = new Blob([data], { type: 'text/plain;charset=utf-8' });
        saveAs(blob, `SQrutiny_${this.construct.name || 'untitled'}.${ext}`);
        this.notify.success(`Exported to ${op}!`);
      } else {
        this.notify.error('Unable to export');
      }
    }
  }

  // TODO:
  // canDeactivate(): Observable<boolean> | boolean {
  // if (this.newExperimentForm.dirty && !this.newExperimentForm.submitted) return confirm('If you leave you\'ll lose all the unsaved data. Are you sure you want to leave this page?'); // Dirty show dialog to user to confirm leaving
  //
  //   return true;
  // }

}
