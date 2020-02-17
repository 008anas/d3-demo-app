import { Component, OnInit, OnDestroy } from '@angular/core';
import { Subscription } from 'rxjs/internal/Subscription';
import { HttpClient, HttpEventType } from '@angular/common/http';
import { FormGroup, FormBuilder } from '@angular/forms';
import { ActivatedRoute, Router } from '@angular/router';
import { finalize } from 'rxjs/operators';

import { NzModalService } from 'ng-zorro-antd/modal';

import { Specie } from '@models/specie';
import { Construct } from '@models/construct';
import { SpecieService } from '@services/specie.service';
import { UserHistory } from 'app/workspace/shared/user-history';

@Component({
  selector: 'sqy-from-file',
  templateUrl: './from-file.component.html',
  styleUrls: ['./from-file.component.scss']
})
export class FromFileComponent implements OnInit, OnDestroy {

  progress = 0;
  form: FormGroup;
  response: any = null;
  sub: Subscription;
  specie: Specie = new Specie();
  species: Specie[] = [];
  isLoading = false;
  history: UserHistory = null;
  isUploading = false;
  construct: Construct = null;
  exampleFile = 'assets/files/NC_044937.gbk';

  constructor(
    private http: HttpClient,
    private route: ActivatedRoute,
    private specieSrvc: SpecieService,
    private formBuilder: FormBuilder,
    private router: Router,
    private modal: NzModalService
  ) { }

  ngOnInit() {
    this.form = this.formBuilder.group({
      file: ['']
    });
    this.sub = this.route.queryParams.subscribe(params => this.specie.slug = params.specie || null);
    this.getSpecie();
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

  getSpecies() {
    this.isLoading = true;
    this.specieSrvc.getAll()
      .pipe(finalize(() => this.isLoading = false))
      .subscribe(
        data =>
          this.species = data.map((e: any) => {
            const specie = new Specie().deserialize(e);
            if (!this.specie.slug && specie.default) {
              this.specie = Object.assign({}, specie);
            }
            return specie;
          })
      );
  }

  onSubmit() {
    const params = new FormData();
    params.append('file', this.form.get('file').value);
    this.response = null;
    this.progress = 0;
    this.isUploading = true;
    this.http.post('http://127.0.0.1:8000/api/v1/optimize_seq/from-file', params, { reportProgress: true, observe: 'events' })
      .pipe(finalize(() => this.isUploading = false))
      .subscribe(
        event => {
          switch (event.type) {
            case HttpEventType.Sent:
              console.log('Request has been made!');
              break;
            case HttpEventType.ResponseHeader:
              console.log('Response header has been received!');
              break;
            case HttpEventType.UploadProgress:
              this.progress = Math.round(event.loaded / event.total * 100);
              console.log(`Uploaded! ${this.progress}%`);
              break;
            case HttpEventType.Response:
              this.construct = new Construct().deserialize(event.body);
              this.modal.success({
                nzTitle: 'Success!',
                nzContent: `${this.history.construct.name} has been submitted: In a few moment you'll be redirected to results page. Please be patient.`
              });
              setTimeout(() => this.router.navigate(['/workspace', this.history.id]), 3000);
          }
        },
        err => {
          this.response = err.msg;
          console.log(err);
        }
      );
  }

  upload(file) {
    console.log(file)
  }

}
