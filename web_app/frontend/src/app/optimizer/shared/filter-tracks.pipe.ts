import { Pipe, PipeTransform } from '@angular/core';
import { Track } from './track';

@Pipe({
  name: 'filterTracks'
})
export class FilterTracksPipe implements PipeTransform {

  transform(items: Track[], query: string): Track[] {
    if (!items) { return []; }
    if (!query) { return items; }
    query = query.toLowerCase();
    return items.filter(item => [item.label || '', item.type || '', item.sequence || ''].some(f => f.toLowerCase().includes(query))
    );
  }

}
