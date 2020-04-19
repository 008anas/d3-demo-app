import { Component, OnInit, OnDestroy } from '@angular/core';
import { Router, ActivatedRoute } from '@angular/router';
import { Subscription } from 'rxjs/internal/Subscription';
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
import { NavService } from '@services/nav.service';
import { FeatureService } from '../shared/feature.service';
import { Feature } from '../shared/feature';

const NAME_IN_SESSION = 'sqy_construct';

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

  private sub: Subscription;
  tracks: Track[] = [];
  track: Track = null;
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
  history: UserHistory = null;
  construct: Construct = new Construct();
  showIndexes = true;
  view = 'general';
  locked = false;
  search: string;
  sketcherLoading = false;
  autoSave = true;
  features: Feature[] = [];
  featuresArray: string[];

  constructor(
    private route: ActivatedRoute,
    private specieSrvc: SpecieService,
    private trackSrvc: TrackService,
    private constructSrvc: ConstructService,
    private sqrutinySrvc: SqrutinyService,
    private featureSrvc: FeatureService,
    private fileSrvc: FileService,
    private router: Router,
    private notify: NzMessageService,
    private navSrvc: NavService
  ) { }

  ngOnInit() {
    this.isLoading = true;
    this.sub = this.route.queryParams.subscribe(params => {
      this.specie.slug = params.specie || null;
      params.construct ? this.getConstruct(params.construct) : this.initConstruct();
    });
    if (this.specie.slug) {
      this.getSpecie();
    }
    this.getTracks();
    this.getSpecies();
    this.getFeatures();
  }

  ngOnDestroy() {
    if (this.sub) {
      this.sub.unsubscribe();
    }
  }

  private initConstruct() {
    if (this.getFromSession()) {
      this.construct = this.getFromSession();
    }
  }

  private isConstructEmpty() {
    for (const key in this.construct) {
      if (this.construct[key] && this.construct[key] !== null && this.construct[key] !== '') {
        return false;
      }
    }
  }

  private getSpecie() {
    this.isLoading = true;
    this.specieSrvc
      .getBySlug(this.specie)
      .pipe(finalize(() => this.isLoading = false))
      .subscribe(data => this.specie.deserialize(data));
  }

  private getConstruct(id: any) {
    if (id) {
      this.sketcherLoading = true;
      this.constructSrvc
        .getById(id)
        .pipe(finalize(() => this.sketcherLoading = false))
        .subscribe(
          (data: Construct) => {
            this.construct = new Construct().deserialize(data);
            this.saveInSession();
          }
        );
    }
  }

  getExampleConstruct() {
    this.sketcherLoading = true;
    this.constructSrvc.getExample()
      .pipe(finalize(() => this.sketcherLoading = false))
      .subscribe(
        (data: Construct[]) => {
          if (data.length > 0) {
            this.construct = new Construct().deserialize(data[0]);
            this.saveInSession();
          } else { this.notify.warning('Unable to load model construct'); }
        },
        err => this.notify.warning(err || 'Unable to load model construct')
      );
  }

  getSpecies() {
    this.isLoading = true;
    this.specieSrvc
      .getAll()
      .pipe(finalize(() => this.isLoading = false))
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
      .pipe(finalize(() => this.isTracksLoading = false))
      .subscribe(data => (this.categories = data));
  }

  getTracks() {
    this.isLoading = true;
    this.trackSrvc
      .getAll()
      .pipe(finalize(() => this.isLoading = false))
      .subscribe(data => {
        this.tracks = data;
        this.new();
      });
  }

  private getFeatures() {
    this.isLoading = true;
    this.featureSrvc
      .getAll()
      .pipe(finalize(() => this.isLoading = false))
      .subscribe(data => {
        this.features = data.map((e: any) => new Feature().deserialize(e));
        this.featuresArray = Object.assign([], this.features.map(f => f.alias));
      }
      );
  }

  featChange(alias: string, isChecked: boolean) {
    isChecked ? this.featuresArray.push(alias) : this.featuresArray = this.featuresArray.filter(f => f !== alias);
  }

  new() {
    this.construct.dna_seq = '';
    if (this.tracks.some(e => e.default) && !this.construct.tracks) {
      this.construct.tracks = [
        Object.assign(
          {},
          this.tracks.find(e => e.default)
        )
      ];
      this.saveInSession();
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
      this.removeFromSession();
    }
  }

  removeTrack(x: Track) {
    const i: number = this.construct.tracks.indexOf(x);
    if (i !== -1) {
      this.construct.tracks.splice(i, 1);
      this.saveInSession();
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
            if (!this.construct.tracks) {
              this.construct.tracks = [];
            }
            this.construct.tracks.push(Object.assign({}, e));
          }
        });
      });
      this.saveInSession();
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
      this.sqrutinySrvc.fromConstruct(this.construct, this.featuresArray.length > 0 ? this.featuresArray : null)
        .pipe(finalize(() => this.sketcherLoading = false))
        .subscribe(
          (data: UserHistory) => {
            this.history = new UserHistory().deserialize(data);
            this.navSrvc.updateBadge();
            this.submitted = true;
            this.saveInSession();
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
      this.track.pos = i;
    }
  }

  addTrack(track: Track) {
    if (track.pos > -1) {
      const length = this.construct.dna_seq.length;
      this.construct.tracks[track.pos] = track;
      this.construct.tracks[track.pos].start = length + 1;
      this.construct.tracks[track.pos].end = length + track.sequence.length;
      this.construct.dna_seq += track.sequence;
      document.getElementById('track' + track.pos).classList.remove('invalid'); // Now is valid
      this.saveInSession();
    }
    this.track = null;
  }

  changeTrack(pos: number) {
    if (pos > -1 && this.construct.tracks[pos]) {
      this.openSidebar(this.construct.tracks[pos], pos);
      this.saveInSession();
    }
  }

  moveTrack(x: number, i: number) {
    const pos = x + i;
    if (-1 < pos && pos <= this.construct.tracks.length - 1) {
      this.construct.tracks = Utils.array_move(this.construct.tracks, x, pos);
      this.saveInSession();
    }
  }

  checkTracks() {
    let flag = true;
    this.construct.tracks.map((t, i) => {
      if (!t.sequence) {
        document.getElementById('track' + i).classList.add('invalid');
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

  // Session storage manage
  private getFromSession() {
    return JSON.parse(sessionStorage.getItem(NAME_IN_SESSION)) || null;
  }

  saveInSession() {
    if (this.autoSave && this.construct) {
      sessionStorage.setItem(NAME_IN_SESSION, JSON.stringify(this.construct));
    }
  }

  private removeFromSession() {
    sessionStorage.removeItem(NAME_IN_SESSION);
  }
}
