import { Component, OnInit, OnDestroy } from '@angular/core';
import { Router, ActivatedRoute } from '@angular/router';
import { Subscription } from 'rxjs/internal/Subscription';
import { Observable } from 'rxjs';
import { finalize } from 'rxjs/operators';

import { NzMessageService } from 'ng-zorro-antd/message';

import { SpecieService } from '@services/specie.service';
import { Specie } from '@models/specie';
import { Track } from '../shared/track';
import { TrackService } from '../shared/track.service';
import { ConstructService } from '@services/construct.service';
import { UserHistory } from 'app/workspace/shared/user-history';
import Utils from 'app/shared/utils';
import { Construct } from '@models/construct';
import { SqrutinyService } from '@services/sqrutiny.service';
import { FileService } from '@services/file.service';

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
  history: UserHistory = null;
  trackHovered: Track = null;
  hoveredName: string = null;
  isLoading = false;
  response: any = null;
  categories: Category[] = [];
  specie: Specie = new Specie();
  species: Specie[] = [];
  submitted = false;
  showPicker = false;
  zoom = 75;
  isTracksLoading = false;
  construct: Construct = new Construct();
  showIndexes = true;
  view = 'general';
  locked = false;
  search: string;
  sketcherLoading = false;

  constructor(
    private route: ActivatedRoute,
    private specieSrvc: SpecieService,
    private trackSrvc: TrackService,
    private constructSrvc: ConstructService,
    private sqrutinySrvc: SqrutinyService,
    private fileSrvc: FileService,
    private router: Router,
    private notify: NzMessageService
  ) {
    this.construct.tracks = [];
    this.construct.dna_seq = '';
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
    if (this.sub) {
      this.sub.unsubscribe();
    }
  }

  getSpecie() {
    this.isLoading = true;
    this.specieSrvc
      .getBySlug(this.specie)
      .pipe(finalize(() => (this.isLoading = false)))
      .subscribe(data => this.specie.deserialize(data));
  }

  getConstruct() {
    if (this.construct.id) {
      this.sketcherLoading = true;
      this.constructSrvc
        .getById(this.construct.id)
        .pipe(finalize(() => (this.sketcherLoading = false)))
        .subscribe(
          data => (this.construct = new Construct().deserialize(data))
        );
    }
  }

  getExampleConstruct() {
    this.sketcherLoading = true;
    this.constructSrvc.getExample()
      .subscribe(
        data =>
          data.length
            ? (this.construct = new Construct().deserialize(data[0]))
            : this.notify.warning('Unable to load model construct'),
        err => this.notify.warning(err || 'Unable to load model construct'),
        () => (this.sketcherLoading = false)
      );
  }

  getSpecies() {
    this.isLoading = true;
    this.specieSrvc
      .getAll()
      .pipe(finalize(() => (this.isLoading = false)))
      .subscribe(
        data =>
          (this.species = data.map((e: any) => {
            const specie = new Specie().deserialize(e);
            if (specie.default && !this.specie.slug) {
              this.specie = Object.assign({}, specie);
            }
            return specie;
          }))
      );
  }

  getTracksByCategories() {
    this.isTracksLoading = true;
    this.trackSrvc
      .getByCategories()
      .pipe(finalize(() => (this.isTracksLoading = false)))
      .subscribe(data => (this.categories = data));
  }

  getTracks() {
    this.isLoading = true;
    this.trackSrvc
      .getAll()
      .pipe(finalize(() => (this.isLoading = false)))
      .subscribe(data => {
        this.tracks = data;
        this.new();
      });
  }

  new() {
    if (this.tracks.some(e => e.default)) {
      this.construct.tracks = [
        Object.assign(
          {},
          this.tracks.find(e => e.default)
        )
      ];
    }
  }

  clear() {
    if (
      confirm(
        'Are you sure you want to clear all tracks from construct? These action cannot be reverted.'
      )
    ) {
      this.construct.tracks = [];
      this.construct.dna_seq = '';
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
      if (
        !confirm(
          'You really want to load an example construct? You\'re gonna lose all actual data. Proceed?'
        )
      ) {
        flag = false;
      }
    }

    if (flag) {
      this.getExampleConstruct();
    }
  }

  submit() {
    if (this.checkTracks()) {
      this.response = null;
      this.sketcherLoading = true;
      this.construct.specie_tax_id = this.specie.tax_id;
      this.sqrutinySrvc.fromSketch(this.construct)
        .pipe(finalize(() => (this.sketcherLoading = false)))
        .subscribe(
          (data: UserHistory) => {
            this.history = new UserHistory().deserialize(data);
            this.submitted = true;
            setTimeout(() => {
              this.submitted = false;
              this.router.navigate(['/workspace', this.history.id]);
            }, 3000);
          },
          err => this.notify.error(err)
        );
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
      this.construct.tracks[track['pos']].start =
        this.construct.dna_seq.length + 1;
      this.construct.tracks[track['pos']].end =
        this.construct.dna_seq.length + track.sequence.length;
      this.construct.dna_seq += track.sequence;
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
    if (-1 < pos && pos <= this.construct.tracks.length - 1) {
      this.construct.tracks = Utils.array_move(this.construct.tracks, x, pos);
    }
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

  toggleSelection(event: { target: { checked: boolean } }) {
    this.categories.map(c =>
      c.elements.map(e => (e.selected = event.target.checked))
    );
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
        case 'Excel':
          // data = this.construct.tracks;
          ext = 'xlsx';
          break;
        case 'JSON':
          data = JSON.stringify({ construct: this.construct });
          ext = 'json';
          break;
      }

      if (data) {
        this.fileSrvc.saveFileAs(data, 'text/plain;charset=utf-8', `SQrutiny_${this.construct.name || 'untitled'}.${ext}`);
        this.notify.success(`Exported to ${op}!`);
      } else {
        this.notify.error('Unable to export');
      }
    }
  }

  // TODO:
  canDeactivate(): Observable<boolean> | boolean {
    // if (!this.submitted && this.construct.tracks.length) return confirm('If you leave you\'ll lose all the unsaved data. Are you sure you want to leave this page?'); // Dirty show dialog to user to confirm leaving

    return true;
  }
}
