import { Component, OnInit, OnDestroy } from '@angular/core';
import { ActivatedRoute } from '@angular/router';
import { Subscription } from 'rxjs/internal/Subscription';
import { finalize } from 'rxjs/operators';

import { Construct } from '../construct/shared/construct';
import { ConstructService } from '../construct/shared/construct.service';

@Component({
  selector: 'sqy-editor',
  templateUrl: './editor.component.html',
  styleUrls: ['./editor.component.scss']
})
export class EditorComponent implements OnInit, OnDestroy {

  sub: Subscription;
  construct: Construct = null;
  isLoading = true;

  constructor(
    private route: ActivatedRoute,
    private constructSrvc: ConstructService
  ) { }

  ngOnInit() {
    this.construct = new Construct().deserialize(this.route.snapshot.data.construct);

    if (this.construct) {
      this.loadConstructInEditor();
    } else {
      this.sub = this.route.queryParams.subscribe(params => this.construct.id = params.construct || null);

      if (this.construct.id) this.getConstruct();
    }
  }

  ngOnDestroy() {
    if (this.sub) this.sub.unsubscribe();
  }

  getConstruct() {
    this.isLoading = true;
    this.constructSrvc.getById(this.construct.id)
      .pipe(finalize(() => this.isLoading = false))
      .subscribe(data => {
        this.construct.deserialize(data);
        this.loadConstructInEditor();
      }
      );
  }

  loadConstructInEditor() {
    const editor = window.createVectorEditor(document.getElementById('vector_editor') || 'createDomNodeForMe'); /* createDomNodeForMe will make a dom node for you and append it to the document.body*/
    editor.updateEditor({
      sequenceData: {
        sequence: this.construct.dna_seq,
        circular: this.construct.circular,
        sequenceFileName: 'pj5_00001.gb',
        name: 'pj5_00001',
        size: 8832,
        description: '',
        translations: {
          '55a4a061f0c5b500012a8qqqq': {
            name: 'Operator I2 and I1',
            type: 'protein_bind',
            id: '55a4a061f0c5b500012a8qqqq',
            start: 1123,
            end: 1162,
            strand: 1,
            notes: [],
            forward: true,
            annotationType: 'translation'
          },
          '55a4a061f0c5b500312340a8qqqq': {
            name: 'pBAD promoter',
            type: 'promoter',
            id: '55a4a061f0c5b500312340a8qqqq',
            start: 1163,
            end: 1188,
            strand: 1,
            notes: [],
            forward: true,
            annotationType: 'translation'
          }
        },
        features: this.construct.tracks.map(t => {
          return {
            color: t.color, //you can override the default color for each individual feature if you want
            name: t.label,
            type: t.type,
            start: t.start, //start and end are 0-based inclusive for all annotations
            end: t.end,
            forward: true //ie true=positive strand     false=negative strange
          }
        })
        // {
        //   '55a4a061f0c5b50asd00a8bfaf5': {
        //     name: 'Operator I2 and I1',
        //     type: 'protein_bind',
        //     id: '55a4a061f0c5b50asd00a8bfaf5',
        //     start: 4,
        //     end: 20,
        //     strand: 1,
        //     notes: [],
        //     color: '#BBBBBB',
        //     forward: true,
        //     annotationType: 'feature'
        //   },
        //   '55a4a061f0c5b5000a8bfaf8': {
        //     name: 'GFPuv',
        //     type: 'CDS',
        //     id: '55a4a061f0c5b5000a8bfaf8',
        //     start: 1235,
        //     end: 2018,
        //     strand: 1,
        //     notes: [
        //       {
        //         name: 'vntifkey',
        //         value: '4',
        //         quoted: true
        //       }
        //     ],
        //     color: '#BBBBBB',
        //     forward: true,
        //     annotationType: 'feature'
        //   },
        // }
        ,
        parts: {},
        primers: {},
        fromFileUpload: false
      },
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
          },
          // ...etc
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
        dnaOrAA: 'DNA', //or 'AA'
        ambiguousOrLiteral: 'LITERAL', //or 'AMBIGUOUS'
        highlightAll: false,
        matchNumber: 0
      },
      deletionLayers: {},
      replacementLayers: {},
      instantiated: true
    });
  }

}
