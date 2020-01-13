import { Component, OnInit } from '@angular/core';
import { ActivatedRoute, Router } from '@angular/router';
import { HttpClient } from '@angular/common/http';

import { environment as env } from 'src/environments/environment';
import { UserHistory } from '../shared/user-history';
import { HistoryService } from '../shared/history.service';
import { JobService } from '../shared/job.service';
import { NotifyService } from 'src/app/shared/services/notify.service';
import { Track } from 'src/app/optimizer/shared/track';
import Utils from 'src/app/shared/utils';
import { Job } from '../shared/job';


@Component({
  selector: 'sqy-history',
  templateUrl: './history.component.html',
  styleUrls: ['./history.component.scss']
})
export class HistoryComponent implements OnInit {

  history: UserHistory = null;
  isLoading = false;
  options: any[] = [];
  isSearching = false;
  msg: string = null;

  constructor(
    private route: ActivatedRoute,
    private router: Router,
    private historySrvc: HistoryService,
    private http: HttpClient,
    private jobSrvc: JobService,
    private notify: NotifyService
  ) { }

  ngOnInit() {
    this.history = new UserHistory().deserialize(this.route.snapshot.data.history);
    this.getJob();
  }

  getJob() {
    this.isLoading = true;
    this.options = [];
    this.jobSrvc.getById(this.history.job_id)
      .subscribe(
        (data: Job) => {
          this.history.job = data;
          // this.history.job.status = 'started';
          // if (this.history.isDone()) {
          //   for (let key in this.history.job.result) {
          //     let value = this.history.job.result[key];
          //     this.options.push({
          //       name: key,
          //       display: true,
          //       data: value
          //     });
          //   };
          //   setTimeout(() => {
          //     document.querySelectorAll('.protvista').forEach((x: any) => x.setAttribute('length', this.history.construct.dna_seq.length));
          //     document.querySelector('#dna_sequence')['data'] = this.history.construct.dna_seq;
          //     document.querySelector('#protein_sequence')['data'] = ' ' + this.history.construct.protein_seq.split('').join('  ') + ' ';
          //     document.querySelector('#interpro-track-residues')['data'] = this.getTrackView(this.history.construct.tracks);
          //     document.querySelectorAll('.var-graph').forEach((x: any, index: number) => x.data = JSON.parse(this.options[index].data));
          //   }, 100);
          // }
        },
        err => {
          console.log(err);
        },
        () => this.isLoading = false
      );
  }

  getStatusClass() {
    return {
      success: this.history.isDone(),
      active: this.history.isActive(),
      error: !this.history.isActive() && !this.history.isDone()
    };
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

  getTrackView(tracks: Track[]) {
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

  clearHighlight() {
    document.querySelectorAll('.protvista').forEach((x: any) => x.fixedHighlight = undefined);
  }

  searchMotif(event: any) {
    if (event.target.value) {
      this.isSearching = true;
      this.msg = null;
      const params = { params: { sequence: this.history.construct.dna_seq, motif: event.target.value } };
      this.http.get(`${env.endpoints.api}/search-motif`, params)
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
              this.msg = 'No match was found';
            }
          },
          err => {
            this.msg = err.msg;
          },
          () => this.isSearching = false
        );
    }
  }

  deleteHistory() {
    if (confirm('Are you sure?')) {
      this.isLoading = true;
      this.historySrvc.delete(this.history.id)
        .subscribe(() => {
          this.notify.success('History deleted!', 'bottom-right', true);
          this.router.navigate(['/workspace']);
        },
          err => this.notify.error(err.msg || 'Unable to delete history', 'bottom-right'),
          () => this.isLoading = false);
    }
  }

}
