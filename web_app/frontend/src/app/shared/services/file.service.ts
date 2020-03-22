import { Injectable } from '@angular/core';

import { saveAs } from 'file-saver';
import * as XLSX from 'xlsx';

import { main } from '@config/main';

const EXCEL_TYPE = 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet;charset=UTF-8';
const EXCEL_EXTENSION = '.xlsx';

class ExcelData {
  name: string;
  data: any;
}

@Injectable({
  providedIn: 'root'
})
export class FileService {

  constructor() { }

  public exportAsExcelFile(data: ExcelData[], fileName: string): void {
    var wb: XLSX.WorkBook = { SheetNames: [], Sheets: {} };
    data.forEach(d => {
      const ws: XLSX.WorkSheet = XLSX.utils.json_to_sheet(d.data);
      wb.SheetNames.push(d.name);
      wb.Sheets[d.name] = ws;
    });
    const excelBuffer: any = XLSX.write(wb, { bookType: 'xlsx', type: 'array' });
    this.saveFileAs(excelBuffer, EXCEL_TYPE, `${main.appName}_${fileName}_export${EXCEL_EXTENSION}`);
  }

  public saveFileAs(buffer: any, type: string, fileName: string): void {
    const data: Blob = new Blob([buffer], { type });
    saveAs(data, fileName);
  }
}
