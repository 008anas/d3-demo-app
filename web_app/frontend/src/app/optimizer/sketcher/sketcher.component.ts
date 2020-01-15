import { Component, OnInit, HostListener, OnDestroy } from '@angular/core';
import { Router, ActivatedRoute } from '@angular/router';
import { Subscription } from 'rxjs/internal/Subscription';
import { HttpClient } from '@angular/common/http';
import { finalize } from 'rxjs/operators';

import { anyToJson, jsonToGenbank, jsonToFasta } from "bio-parsers";
import FileSaver from "file-saver";

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
  isLoading = false;
  response: any = null;
  categories: Category[] = [];
  specie: Specie = new Specie();
  species: Specie[] = [];
  isSubmitting = false;
  submitted = false;
  showPicker = false;
  zoom: number = 75;
  isTracksLoading = false;
  construct: Construct = new Construct();
  showIndexes = true;
  totalSeq = '';

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
  }

  ngOnInit() {
    this.isLoading = true;
    this.sub = this.route.queryParams.subscribe(params => {
      this.construct.id = params.uuid || null;
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
    this.sub.unsubscribe();
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
      .subscribe(data => this.construct = new Construct().deserialize(data));
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
      this.totalSeq = '';
    }
  }

  removeElem(x: Track) {
    const i: number = this.construct.tracks.indexOf(x);
    if (i !== -1) {
      this.construct.tracks.splice(i, 1);
    }
  }

  someSelected() {
    return this.categories.some(c => c.elements.some(t => t.selected));
  }

  addTracks() {
    if (this.someSelected()) {
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
    this.construct.label = 'Example construct';
    this.construct.tracks = Object.assign([], this.tracks.slice(0, 4));
    this.construct.tracks[0].label = 'First element';
    this.construct.tracks[0].sequence = 'CGGCGAGCGGCGAGTAACGGCGAGCGGCGAGTAAAATATATAAAATGAGCGGAGAGCGCGCGGCGAGCGGCGAGTAACGGCGAGCGGCGAGTAAAATATATAAAATGAGCGGAGAGCGCG';
    this.construct.tracks[0].color = '#96c582';
    this.totalSeq += this.construct.tracks[0].sequence;
    this.construct.tracks[1].label = 'Second element';
    this.construct.tracks[1].sequence = 'CGGCGAGCGGCGAGTAACGGCGAGCGGCGAGTAAAATATATAAAATGAGCGGAGAGCGCGCGGCGAGCGGCGAGTAACGGCGAGCGGCGAGTAAAATATATAAAATGAGCGGAGAGCGCG';
    this.construct.tracks[1].color = '#cf0500';
    this.totalSeq += this.construct.tracks[1].sequence;
    this.construct.tracks[2].label = 'Third element';
    this.construct.tracks[2].sequence = 'CGGCGAGCGGCGAGTAACGGCGAGCGGCGAGTAAAATATATAAAATGAGCGGAGAGCGCGCGGCGAGCGGCGAGTAACGGCGAGCGGCGAGTAAAATATATAAAATGAGCGGAGAGCGCG';
    this.construct.tracks[2].color = '#4b4fce';
    this.totalSeq += this.construct.tracks[2].sequence;
    this.construct.tracks[3].label = 'Fourth element';
    this.construct.tracks[3].sequence = 'CGGCGAGCGGCGAGTAACGGCGAGCGGCGAGTAAAATATATAAAATGAGCGGAGAGCGCGCGGCGAGCGGCGAGTAACGGCGAGCGGCGAGTAAAATATATAAAATGAGCGGAGAGCGCG';
    this.construct.tracks[3].color = '#f0ca68';
    this.totalSeq += this.construct.tracks[3].sequence;
  }

  submit() {
    if (this.checkTracks()) {
      this.response = null;
      this.isSubmitting = true;
      this.construct['specie_tax_id'] = this.specie.tax_id;
      this.construct.tracks.map(t => t['element'] = t.id);
      this.http.post(`${env.endpoints.api}/optimize_seq/from-sketch`, this.construct)
        .pipe(finalize(() => this.isSubmitting = false))
        .subscribe(
          (data: UserHistory) => {
            this.history = new UserHistory().deserialize(data);
            this.submitted = true;
          },
          err => {
            this.notify.error(err.msg, 'bottom-right');
          });
    }
  }

  // Track Details Sidebar

  openSidebar(e: Track, i: number) {
    this.track = Object.assign({}, e);
    this.track['pos'] = i;
  }

  addTrack(track: Track) {
    if (track['pos'] > -1) {
      this.construct.tracks[track['pos']] = track;
      this.totalSeq += track.sequence;
      delete this.construct.tracks[track['pos']]['invalid']; // Now is valid
    }
    this.track = null;
  }

  changeTrack(pos: number) {
    if (pos > -1 && this.construct.tracks[pos]) {
      this.openSidebar(this.construct.tracks[pos], pos);
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

  toggleSelection(event: { target: { checked: boolean; }; }) {
    this.categories.map(c => c.elements.map(e => e.selected = event.target.checked));
  }

  // Export / Save Construct

  downloadAs(op: string) {
    if (this.construct.tracks.length) {

      let data, ext;

      switch (op) {
        case 'gb':
          data = jsonToGenbank(this.construct);
          ext = 'gbk';
          break;
        case 'fasta':
          data = jsonToFasta(this.construct);
          ext = 'fasta'
          break;
        case 'json':
          data = JSON.stringify({ 'construct': this.construct });
          ext = 'json';
          break;
      }

      const blob = new Blob([data], { type: "text/plain" });
      const filename = `SQrutiny_${this.construct.label || "untitled"}.${ext}`;
      FileSaver.saveAs(blob, filename);
    }
  }



  goToHistory() {
    if (this.history.id) {
      this.router.navigate(['/workspace', this.history.id]);
      this.submitted = false;
    }
  }

  // TODO:
  // canDeactivate(): Observable<boolean> | boolean {
  // if (this.newExperimentForm.dirty && !this.newExperimentForm.submitted) return confirm('If you leave you\'ll lose all the unsaved form_data. Are you sure you want to leave this page?'); // Dirty show dialog to user to confirm leaving
  //
  //   return true;
  // }

}
