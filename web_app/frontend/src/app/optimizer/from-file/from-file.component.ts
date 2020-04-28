import { Component, OnInit, OnDestroy } from '@angular/core';
import { Subscription } from 'rxjs/internal/Subscription';
import { ActivatedRoute, Router } from '@angular/router';
import { HttpClient, HttpResponse, HttpEvent, HttpEventType, HttpRequest } from '@angular/common/http';
import { finalize } from 'rxjs/operators';

import { UploadXHRArgs } from 'ng-zorro-antd/upload';
import { NzModalService } from 'ng-zorro-antd/modal';
import { NzMessageService } from 'ng-zorro-antd/message';

import { environment as env } from '@env/environment';
import { Specie } from '@models/specie';
import { Construct } from '@models/construct';
import { SpecieService } from '@services/specie.service';
import { UserHistory } from 'app/workspace/shared/user-history';
import { TextModalComponent } from '@components/text-modal/text-modal.component';
import { SqrutinyService } from '@services/sqrutiny.service';
import { LoaderService } from '@services/loader.service';
import { NavService } from '@services/nav.service';
import { FeatureService } from '../shared/feature.service';
import { Feature } from '../shared/feature';
import { TrackDisplayComponent } from '@components/track-display/track-display.component';
import { Track } from '../shared/track';

@Component({
  selector: 'sqy-from-file',
  templateUrl: './from-file.component.html',
  styleUrls: ['./from-file.component.scss']
})
export class FromFileComponent implements OnInit, OnDestroy {

  response: any = null;
  sub: Subscription;
  specie: Specie = new Specie();
  species: Specie[] = [];
  isLoading = false;
  construct: Construct = null;
  history: UserHistory = null;
  fileList = [];
  endpoint: string;
  view = true;
  search: string;
  isSubmitting = false;
  allowedExt = ['.gb', '.gbk', '.genbank'];
  isFeaturesLoading = false;
  features: Feature[] = [];
  featuresArray: string[];
  tplModal: any;

  customReq = (item: UploadXHRArgs) => {
    this.loader.startLoading();
    const formData = new FormData();
    formData.append('file', item.file as any);
    const req = new HttpRequest('POST', item.action!, formData, {
      reportProgress: true,
      withCredentials: true
    });
    this.fileList = [item.file];
    this.response = null;
    this.construct = null;
    return this.http.request(req)
      .pipe(finalize(() => this.loader.stopLoading()))
      .subscribe(
        (event: HttpEvent<any>) => {
          if (event.type === HttpEventType.UploadProgress) {
            if (event.total! > 0) {
              (event as any).percent = (event.loaded / event.total!) * 100;
            }
            item.onProgress!(event, item.file!);
          } else if (event instanceof HttpResponse) {
            item.file.status = 'done';
            item.onSuccess!(event.body, item.file!, event);
            this.construct = new Construct().deserialize(event.body);
            this.notify.success('GenBank file parsed successfully!');
          }
        },
        err => this.response = err
      );
  }

  constructor(
    private route: ActivatedRoute,
    private specieSrvc: SpecieService,
    private sqrutinySrvc: SqrutinyService,
    private router: Router,
    private notify: NzMessageService,
    private http: HttpClient,
    private modal: NzModalService,
    private loader: LoaderService,
    private navSrvc: NavService,
    private featureSrvc: FeatureService
  ) {
    this.endpoint = env.endpoints.api + '/constructs/from-genbank';
  }

  ngOnInit() {
    this.sub = this.route.queryParams.subscribe(params => {
      if (params.specie) {
        this.specie.slug = params.specie || null;
        this.getSpecie();
      }
    });
    this.getSpecies();
    this.getFeatures();
  }

  ngOnDestroy() {
    if (this.sub) {
      this.sub.unsubscribe();
    }
    if (this.tplModal) {
      this.tplModal.destroy();
    }
  }

  private getSpecie() {
    this.isLoading = true;
    this.specieSrvc
      .getBySlug(this.specie)
      .pipe(finalize(() => this.isLoading = false))
      .subscribe(data => {
        this.specie.deserialize(data);
        this.getFeatures();
      });
  }

  getSpecies() {
    this.isLoading = true;
    this.specieSrvc
      .getAll()
      .pipe(finalize(() => this.isLoading = false))
      .subscribe(
        data => {
          this.species = data.map((e: any) => {
            const specie = new Specie().deserialize(e);
            if (specie.default && !this.specie.slug) {
              this.specie = Object.assign({}, specie);
            }
            return specie;
          });
          this.getFeatures();
        }
      );
  }

  getFeatures() {
    this.isFeaturesLoading = true;
    this.featureSrvc
      .getAll(this.specie.tax_id || null)
      .pipe(finalize(() => this.isFeaturesLoading = false))
      .subscribe(data => {
        this.features = data.map((e: any) => new Feature().deserialize(e));
        this.featuresArray = Object.assign([], this.features.map(f => f.alias));
      });
  }

  submit() {
    this.response = null;
    this.isSubmitting = true;
    this.construct.specie_tax_id = this.specie.tax_id;
    this.sqrutinySrvc.fromConstruct(this.construct, this.featuresArray.length > 0 ? this.featuresArray : null)
      .subscribe(
        (data: UserHistory) => {
          this.navSrvc.updateBadge();
          this.notify.loading('Your job has been submitted! In a moment you will be redirected. Please be patient');
          setTimeout(() => {
            this.router.navigate(['/workspace', data.id]);
            this.isSubmitting = false;
          }, 3000);
        },
        err => {
          this.isSubmitting = false;
          this.notify.error(err);
        }
      );
  }

  showTrack(track: Track) {
    if (track) {
      this.tplModal = this.modal.create({
        nzContent: TrackDisplayComponent,
        nzComponentParams: { track },
        nzCancelText: null,
        nzOnOk: () => {
          if (this.tplModal) {
            this.tplModal.destroy();
          }
        }
      });
    }
  }

  featChange(alias: string, isChecked: boolean) {
    isChecked ? this.featuresArray.push(alias) : this.featuresArray = this.featuresArray.filter(f => f !== alias);
  }

  textModal(str: string) {
    this.modal.create({
      nzContent: TextModalComponent,
      nzComponentParams: {
        txt: str
      },
      nzFooter: null
    });
  }

}
