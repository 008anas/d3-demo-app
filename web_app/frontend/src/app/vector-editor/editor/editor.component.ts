import { Component, OnInit, OnDestroy } from '@angular/core';
import { ActivatedRoute } from '@angular/router';
import { Subscription } from 'rxjs/internal/Subscription';
import { finalize } from 'rxjs/operators';

import { NzModalService } from 'ng-zorro-antd/modal';
declare var createVectorEditor: any;

import { TitleService } from '@services/title.service';
import { UserHistory } from 'app/workspace/shared/user-history';
import { HistoryPickerComponent } from '../shared/history-picker/history-picker.component';
import { HistoryService } from 'app/workspace/shared/history.service';

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
  isLoading = false;
  historyModal = null;

  constructor(
    private route: ActivatedRoute,
    private titleSrvc: TitleService,
    private modal: NzModalService,
    private historySrvc: HistoryService
  ) { }

  ngOnInit(): void {
    this.isLoading = true;
    this.loadEditor();
    if (this.route.snapshot.data.history) {
      this.history = new UserHistory().deserialize(this.route.snapshot.data.history);
      this.titleSrvc.setTitle(this.history.name);
      this.loadConstructInEditor();
    } else {
      this.getHistoriesCount();
    }
  }

  ngOnDestroy() {
    if (this.sub) { this.sub.unsubscribe(); }
    if (this.historyModal) { this.historyModal.close(); }
  }

  getHistoriesCount() {
    this.historySrvc.getCount()
      .pipe(finalize(() => this.isLoading = false))
      .subscribe(data => {
        if (Number.parseInt(data.count || 0, 10) > 0) {
          this.historyModal = this.modal.create({
            nzTitle: 'Choose a history to load in the editor',
            nzContent: HistoryPickerComponent,
            nzCancelText: 'Just go to the editor',
            nzOkText: 'Load',
            nzOnOk: (picker: HistoryPickerComponent) => {
              if (picker.selected) {
                this.history = Object.assign({}, picker.selected);
                this.loadConstructInEditor();
              }
            },
            nzClassName: 'center-modal-scroll'
          });
        }
      });
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
        features: this.history.construct.tracks.map(t => ({
          color: t.color,
          name: t.label,
          type: t.type,
          start: t.start,
          end: t.end,
          forward: true // ie true=positive strand     false=negative strange
        })),
        fromFileUpload: false
      }
    });
  }

}
