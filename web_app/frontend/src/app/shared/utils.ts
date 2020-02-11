const nameUtils = {
	/**
	 * Reformat name to replaces whitespace with underscores.
	 * @returns {String} New name.
	 */
  reformatName: (pName: string): string => pName.toString().replace(/ /g, '_'),
};
const StringUtil = {
  /** Trims white space at beginning and end of string
   * @returns {String} line
   */
  trim: (line: string): string => line.replace(/^\s+|\s+$/g, ''),

  /** Trims white space at beginning string
   * @returns {String} line
   */
  ltrim: (line: string): string => line.replace(/^\s+/, ''),

  /** Trims white space at end of string
   * @returns {String} line
   */
  rtrim: (line: string): string => line.replace(/\s+$/, ''),

  /** Pads white space at beginning of string
   * @returns {String} line
   */
  lpad: (line: string, padString, length): string => {
    let str = line;
    while (str.length < length) str = padString + str;
    return str;
  },

  /** Pads white space at end of string
   * @param {String} line
   * @returns {String} line
   */
  rpad: (line: string, padString, length): string => {
    let str = line;
    while (str.length < length) str = str + padString;
    return str;
  }
};

export default class Utils {

  static dnaSeqRegex = /^[CAGTNU\n]+$/i;
  static porteinSeqRegex = /^[ARNDCQEGHILKMFPSTWYVX*]+$/i;
  static urlRegex = /^(?:(?:https?|ftp):\/\/)?(?:(?!(?:10|127)(?:\.\d{1,3}){3})(?!(?:169\.254|192\.168)(?:\.\d{1,3}){2})(?!172\.(?:1[6-9]|2\d|3[0-1])(?:\.\d{1,3}){2})(?:[1-9]\d?|1\d\d|2[01]\d|22[0-3])(?:\.(?:1?\d{1,2}|2[0-4]\d|25[0-5])){2}(?:\.(?:[1-9]\d?|1\d\d|2[0-4]\d|25[0-4]))|(?:(?:[a-z\u00a1-\uffff0-9]-*)*[a-z\u00a1-\uffff0-9]+)(?:\.(?:[a-z\u00a1-\uffff0-9]-*)*[a-z\u00a1-\uffff0-9]+)*(?:\.(?:[a-z\u00a1-\uffff]{2,})))(?::\d{2,5})?(?:\/\S*)?$/;

  static array_move(arr: any[], old_index: number, new_index: number) {
    if (new_index >= arr.length) {
      let k = new_index - arr.length + 1;
      while (k--) {
        arr.push(undefined);
      }
    }
    arr.splice(new_index, 0, arr.splice(old_index, 1)[0]);
    return arr;
  }

  static jsonToGenbank(serSeq, options?) {
    options = options || {};
    options.reformatSeqName = options.reformatSeqName !== false;

    if (!serSeq) { return false; }

    try {
      if (serSeq.isProtein || serSeq.type === 'protein' || serSeq.type === 'AA') {
        serSeq.isProtein = true;
        serSeq.sequence = serSeq.proteinSequence || serSeq.sequence;
        options.isProtein = true;
      }
      const content = null;
      const cutUp = typeof serSeq.sequence === 'string' ? Utils.cutUpStr : Utils.cutUpArray;
      if (!serSeq.sequence) { serSeq.sequence = ''; }

      let lines = [];
      lines.push(Utils.createGenbankLocus(serSeq, options));
      if (serSeq.definition || serSeq.description) {
        lines.push('DEFINITION  ' + (serSeq.definition || serSeq.description));
      }

      if (serSeq.extraLines) {
        lines = lines.concat(serSeq.extraLines);
      }
      if (serSeq.comments) {
        serSeq.comments.forEach((comment: string) => lines.push('COMMENT             ' + comment));
      }
      if (serSeq.teselagen_unique_id) {
        lines.push(
          'COMMENT             teselagen_unique_id: ' + serSeq.teselagen_unique_id
        );
      }
      if (serSeq.library) {
        lines.push('COMMENT             library: ' + serSeq.library);
      }

      return content;

      // serSeq.features = serSeq.features.map().concat(
      //   flatMap(pragmasAndTypes, ({ pragma, type }) => {
      //     return flatMap(serSeq[type], ann => {
      //       if (!isObject(ann)) {
      //         return [];
      //       }
      //       ann.notes = pragma
      //         ? {
      //           ...ann.notes,
      //           pragma: [pragma]
      //         }
      //         : ann.notes;
      //       return ann;
      //     });
      //   })
      // );
      //
      // let printedFeatureHeader;
      // serSeq.features.foreach((feat, index) => {
      //   if (!printedFeatureHeader) {
      //     printedFeatureHeader = true;
      //     lines.push('FEATURES             Location/Qualifiers');
      //   }
      //   lines.push(Utils.featureToGenbankString(feat, options));
      // });
      //
      // lines.push('ORIGIN      ');
      // for (let i = 0; i < serSeq.sequence.length; i = i + 60) {
      //   const line = [];
      //   const ind = StringUtil.lpad('' + (i + 1), ' ', 9);
      //   line.push(ind);
      //
      //   for (let j = i; j < i + 60; j = j + 10) {
      //     // line.push(serSeq.sequence.slice(j,j+10).join(''));
      //     line.push(cutUp(serSeq.sequence, j, j + 10));
      //   }
      //   lines.push(line.join(' '));
      // }
      //
      // lines.push('//');
      //
      // content = lines.join('\r\n');
      // // return cb(err, content);
      // return content;
    } catch (e) {
      console.warn('Error processing sequence << Check jsonToGenbank.js');
      console.warn(serSeq);
      console.warn(e.stack);
      return false;
    }
  }

  static cutUpArray(val, start, end) {
    return val.slice(start, end).join('');
  }

  static cutUpStr(val, start, end) {
    return val.slice(start, end);
  }

  static createGenbankLocus(serSeq, options) {
    if (serSeq.sequence.symbols) {
      serSeq.sequence = serSeq.sequence.symbols.split('');
    }

    let tmp;
    let dnaType;
    if (serSeq.isProtein) {
      dnaType = '';
    } else if (serSeq.type === 'RNA') {
      dnaType = 'RNA';
    } else {
      dnaType = 'DNA';
    }
    const date = Utils.getCurrentDateString();

    let line = StringUtil.rpad('LOCUS', ' ', 12);
    let nameToUse = serSeq.name || 'Untitled_Sequence';
    nameToUse = options.reformatSeqName
      ? nameUtils.reformatName(nameToUse)
      : nameToUse;
    line += StringUtil.rpad(nameToUse, ' ', 16);
    line += ' '; // T.H line 2778 of GenbankFormat.as col 29 space
    line += StringUtil.lpad(String(serSeq.sequence.length), ' ', 11);
    line += serSeq.isProtein ? ' aa ' : ' bp '; // col 41
    // if (strandType !== '') {
    // 	tmp =  strandType + '-';
    // } else {
    tmp = '';
    // }
    line += StringUtil.lpad(tmp, ' ', 3);
    line += StringUtil.rpad(dnaType, ' ', 6);
    line += '  ';

    if (!serSeq.circular || serSeq.circular === '0') {
      line += 'linear  ';
      // line += '        ';
    } else {
      line += 'circular';
    }

    line += ' '; // col 64
    if (serSeq.gbDivision) {
      line += StringUtil.rpad(serSeq.gbDivision || 'SYN', ' ', 10);
    }
    // }
    line += ' '; // col 68
    // DOES NOT PARSE DATE USEFULLY ORIGINALLY!
    line += date;
    // line += '\n';

    return line;
  }

  static featureToGenbankString(feat, options) {
    const lines = [];

    const line = '     ' + StringUtil.rpad(feat.type || 'misc_feature', ' ', 16);
    let locStr = '';

    // for(var i=0;i<feat.locations.length;i++) {
    //	var loc = feat.locations[i];
    //	locStr.push((loc.start+1) + '..' + loc.end);
    // }

    if (feat.locations && feat.locations.length > 1) {
      feat.locations.forEach((loc, i) => {
        locStr +=
          Utils.getProteinStart(
            parseInt(loc.start, 10) + (options.inclusive1BasedStart ? 0 : 1),
            options.isProtein
          ) +
          '..' +
          Utils.getProteinEnd(
            parseInt(loc.end, 10) + (options.inclusive1BasedEnd ? 0 : 1),
            options.isProtein
          );

        if (i !== feat.locations.length - 1) {
          locStr += ',';
        }
      });
      locStr = 'join(' + locStr + ')';
    } else {
      locStr +=
        Utils.getProteinStart(
          parseInt(feat.start, 10) + (options.inclusive1BasedStart ? 0 : 1),
          options.isProtein
        ) +
        '..' +
        Utils.getProteinEnd(
          parseInt(feat.end, 10) + (options.inclusive1BasedEnd ? 0 : 1),
          options.isProtein
        );
    }

    // locStr = locStr.join(',');

    if (feat.strand === -1) {
      locStr = 'complement(' + locStr + ')';
    }

    lines.push(line + locStr);

    lines.push(
      Utils.featureNoteInDataToGenbankString('label', feat.name || 'Untitled Feature')
    );

    let notes = feat.notes;
    if (notes) {
      try {
        if (typeof notes === 'string') {
          try {
            notes = JSON.parse(notes);
          } catch (e) {
            console.warn('Warning: Note incorrectly sent as a string.');
            notes = {}; // set the notes to a blank object
          }
        }
        Object.keys(notes).forEach((key) => {
          if (key === 'color' || key === 'labelColor') { return; } //we'll handle this below
          if (notes[key] instanceof Array) {
            notes[key].forEach((value) => lines.push(Utils.featureNoteInDataToGenbankString(key, value)));
          } else {
            console.warn('Warning: Note object expected array values');
            console.warn(notes);
          }
        });
      } catch (e) {
        console.warn('Warning: Note cannot be processed');
      }
    }
    feat.color = (feat.notes && feat.notes.color) || feat.color;
    feat.labelColor = (feat.notes && feat.notes.labelColor) || feat.labelColor;

    // if (
    //   feat.color &&
    //   color.rgb(feat.color).string() !==
    //   color.rgb(featureColors[feat.type]).string() //don't save a color note if the color is already the same as our defaults
    // ) {
    //   lines.push(Utils.featureNoteInDataToGenbankString('color', feat.color));
    // }
    if (feat.labelColor) {
      lines.push(Utils.featureNoteInDataToGenbankString('labelColor', feat.labelColor));
    }

    return lines.join('\r\n');
  }

  static getCurrentDateString() {
    const date = new Date();
    // date = date.toString().split(' ');
    const day = date[2];
    const month = date[1].toUpperCase();
    const year = date[3];
    return day + '-' + month + '-' + year;
  }

  static getProteinStart(val, isProtein) {
    if (!isProtein) { return val; }
    return Math.floor((val + 2) / 3);
  }
  static getProteinEnd(val, isProtein) {
    if (!isProtein) { return val; }
    return Math.floor(val / 3);
  }

  static featureNoteInDataToGenbankString(name, value) {
    return StringUtil.lpad('/', ' ', 22) + name + '="' + value + '"';
  }

  static jsonToFasta(data: any[]) {
    if (!data) { return null; }

    let fastaString = '';
    data.forEach((c: any) => {
      if (c.sequence) {
        console.log(c.sequence);
        fastaString += `>${c.name || 'Untitled Sequence'}`;
        fastaString += `|${c.sequence.length}`;
        fastaString += c.description ? '|' + c.description : '';
        fastaString += '|' + (c.circular ? 'circular' : 'linear') || 'linear';
        fastaString += '\n';
        fastaString += (c.sequence.match(/.{1,80}/g) || []).join('\n');
      }
    });
    return fastaString;
  }

  getPosFromSeq(main_seq: string, seq: string) {
    if (!main_seq || seq || main_seq.length < seq.length) { return null; }
    return main_seq.indexOf(seq);
  }

}
