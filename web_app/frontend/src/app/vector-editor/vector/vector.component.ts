import { Component, OnInit, OnDestroy } from '@angular/core';
import { ActivatedRoute } from '@angular/router';
import { Subscription } from 'rxjs/internal/Subscription';
import { finalize } from 'rxjs/operators';

declare var createVectorEditor: any;

import { UserHistory } from 'app/workspace/shared/user-history';
import { HistoryService } from 'app/workspace/shared/history.service';
import { TitleService } from '@services/title.service';
import { LoaderService } from '@services/loader.service';
import Utils from 'app/shared/utils';

@Component({
  selector: 'sqy-vector',
  templateUrl: './vector.component.html',
  styleUrls: ['./vector.component.scss']
})
export class VectorComponent implements OnInit, OnDestroy {

  sub: Subscription;
  history: UserHistory = null;
  historyId: string;
  editor: any = null;
  response: any = null;
  uuidRegex = Utils.uuidRegex;

  constructor(
    private route: ActivatedRoute,
    private historySrvc: HistoryService,
    private titleSrvc: TitleService,
    private loader: LoaderService
  ) { }

  ngOnInit() {
    this.loadEditor();
    if (this.route.snapshot.data.history) {
      this.history = new UserHistory().deserialize(this.route.snapshot.data.history);
      this.historyId = this.history.id;
      this.titleSrvc.setTitle(this.history.name);
      this.loadConstructInEditor();
    }
  }

  ngOnDestroy() {
    if (this.sub) { this.sub.unsubscribe(); }
  }

  getHistory() {
    this.loader.startLoading();
    this.response = null;
    this.historySrvc.getByIdNot404(this.historyId)
      .pipe(finalize(() => this.loader.stopLoading()))
      .subscribe(data => {
        this.history = new UserHistory().deserialize(data);
        this.historyId = this.history.id;
        this.loadConstructInEditor();
      },
        err => this.response = err);
  }

  loadEditor() {
    this.editor = createVectorEditor(document.getElementById('vector_editor') || 'createDomNodeForMe'); /* createDomNodeForMe will make a dom node for you and append it to the document.body*/
    this.editor.updateEditor({
      annotationLabelVisibility: {
        features: true,
        parts: true,
        cutsites: true
      },
      annotationVisibility: {
        features: true,
        translations: true,
        parts: true,
        orfs: false,
        orfTranslations: false,
        cdsFeatureTranslations: true,
        axis: true,
        cutsites: true,
        reverseSequence: true
      },
      annotationsToSupport: {
        features: true,
        translations: true,
        parts: true,
        orfs: true,
        cutsites: true,
        primers: true,
      },
      panelsShown: [
        [
          {
            id: 'sequence',
            name: 'Sequence Map',
            active: true
          }
        ],
        [
          {
            id: 'circular',
            name: 'Circular',
            active: true
          },
          {
            id: 'rail',
            name: 'Linear',
            active: false
          },
          {
            id: 'properties',
            name: 'Properties',
            active: false
          }
        ]
      ],
      restrictionEnzymes: {
        filteredRestrictionEnzymes: [
          {
            value: 'single',
            label: 'Single cutters',
            cutsThisManyTimes: 1
          }
        ],
        allRestrictionEnzymes: {
          aatii: {
            name: 'AatII',
            site: 'gacgtc',
            forwardRegex: 'gacgtc',
            reverseRegex: 'gacgtc',
            topSnipOffset: 5,
            bottomSnipOffset: 1,
            usForward: 0,
            usReverse: 0,
            color: '#059369'
          },
          acci: {
            name: 'AccI',
            site: 'gtmkac',
            forwardRegex: 'gt[acm][gkt]ac',
            reverseRegex: 'gt[acm][gkt]ac',
            topSnipOffset: 2,
            bottomSnipOffset: 4,
            usForward: 0,
            usReverse: 0,
            color: '#0d994a'
          }
        }
      },
      selectedAnnotations: {
        idMap: {},
        idStack: []
      },
      minimumOrfSize: 300,
      hoveredAnnotation: '',
      caretPosition: -1,
      selectionLayer: {
        start: -1,
        end: -1
      },
      readOnly: false,
      findTool: {
        isOpen: false,
        searchText: '',
        dnaOrAA: 'DNA', // or 'AA'
        ambiguousOrLiteral: 'LITERAL', // or 'AMBIGUOUS'
        highlightAll: false,
        matchNumber: 0
      },
      deletionLayers: {},
      replacementLayers: {},
      instantiated: true
    });
  }

  loadConstructInEditor() {
    this.editor.updateEditor({
      sequenceData: {
        sequence: this.history.construct.dna_seq,
        circular: this.history.construct.circular,
        sequenceFileName: 'pj5_00001.gb',
        size: 8832,
        description: this.history.construct.description || '',
        features: this.history.construct.tracks.map(t => {
          return {
            color: t.color,
            name: t.label,
            type: t.type,
            start: t.start,
            end: t.end,
            forward: true // ie true=positive strand     false=negative strange
          };
        }),
        fromFileUpload: false
      }
    });
  }
}
