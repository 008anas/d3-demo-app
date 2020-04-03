import { Component, OnInit, OnDestroy } from '@angular/core';
import { ActivatedRoute, Router } from '@angular/router';
import { Subscription } from 'rxjs/internal/Subscription';

import { NzModalService } from 'ng-zorro-antd/modal';

import { TitleService } from '@services/title.service';
import { UserHistory } from 'app/workspace/shared/user-history';
import { HistoryPickerComponent } from '../shared/history-picker/history-picker.component';

declare var createVectorEditor: any;

@Component({
  selector: 'sqy-editor',
  templateUrl: './editor.component.html',
  styleUrls: ['./editor.component.scss']
})
export class EditorComponent implements OnInit, OnDestroy {

  sub: Subscription;
  history: UserHistory = null;
  editor: any = null;
  response: any = null;

  constructor(
    private route: ActivatedRoute,
    private titleSrvc: TitleService,
    private modal: NzModalService,
    private router: Router
  ) {
    this.router.routeReuseStrategy.shouldReuseRoute = () => false;
  }

  ngOnInit(): void {
    this.loadEditor();
    if (this.route.snapshot.data.history) {
      this.history = new UserHistory().deserialize(this.route.snapshot.data.history);
      this.titleSrvc.setTitle(this.history.name);
      this.loadConstructInEditor();
    } else {
      this.modal.create({
        nzTitle: 'Choose a history to load in the editor',
        nzContent: HistoryPickerComponent,
        nzCancelText: 'Go to editor without loading construct',
        nzOkText: 'Load',
        nzOnOk: (data: HistoryPickerComponent) => {
          if (data.selected) {
            this.history = data.selected;
            this.loadConstructInEditor();
          }else{
            this.loadEditor();
          }
        },
        nzClassName: 'center-modal-scroll'
      });
    }
  }


  ngOnDestroy() {
    if (this.sub) { this.sub.unsubscribe(); }
  }

  loadEditor() {
    this.editor = createVectorEditor(document.getElementById('vector_editor') || 'createDomNodeForMe');
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
        allRestrictionEnzymes: {}
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
