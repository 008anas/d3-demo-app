import ProtvistaTrack from "protvista-track";
import {
  scaleLinear,
  line,
  curveMonotoneX,
  select,
  max,
  extent
} from "d3";

class SqrutinyMatrix extends ProtvistaTrack {

  constructor() {
    super();
    this._line = line()
      .defined(d => !isNaN(d.y))
      .x(d => this.getXFromSeqPosition(d.x))
      .y(d => this._yScale(d.y))
      .curve(curveMonotoneX);
  }

  connectedCallback() {
    super.connectedCallback();

    this._height = parseInt(this.getAttribute("height")) || 40;
    this._color = this.getAttribute("color") || 'darkgrey';
    this._cutoffs = this.getAttribute("cutoffs") || [];
    this._yScale = scaleLinear();
    this._xExtent;
    this._yExtent;
  }

  set data(data) {
    this._data = undefined;
    if (!Array.isArray(data) || data.length < 1) {
      return;
    }
    this._data = data.map(d => {
      return {
        x: d.pos,
        y: d.score
      };
    });
    this._createTrack();
  }

  set color(color) {
    this._color = color;
  }

  set cutoffs(cutoffs) {
    if (!Array.isArray(cutoffs) || cutoffs.length < 1) {
      this._cutoffs = [];
    }
    this._cutoffs = cutoffs;
    this.refresh();
  }

  updateYScale() {
    this._xExtent = extent(this._data, d => parseInt(d.x));
    this._yExtent = extent(this._data, d => d.y);

    this.xScale.domain(this._xExtent).range([0, this.width]);
    this._yScale.domain(this._yExtent).range([this._height - this.margin.bottom, this.margin.top]);
  }

  _createTrack() {
    if (!this.svg) {
      this.svg = select(this)
        .append('div')
        .append('svg')
        .attr('width', this.getWidthWithMargins())
        .attr('height', this._height);

      this.cutoff_g = this.svg.append('g').attr('class', 'cutoff');

      this.trackHighlighter.appendHighlightTo(this.svg);
    }

    // Create the visualisation here
    this.updateYScale();
    this.refresh();
  }

  refresh() {
    if (!this.svg) return;
    this.svg.selectAll('path').remove();
    this.svg.selectAll('line').remove();
    this.svg.append('path')
      .datum(this._data)
      .attr('d', this._line)
      .attr('fill', 'none')
      .attr('stroke', this._color)
      .attr('stroke-width', `1.2px`);
    if (this._cutoffs && this._cutoffs.length) {
      this.cutoff_g
        .selectAll('g')
        .data(this._cutoffs)
        .enter()
        .append('line')
        .style('stroke', 'rgb(0, 0, 255)')
        .style('stroke-dasharray', '5px')
        .style('stroke-width', '1.5px')
        .attr('x1', d => this.getXFromSeqPosition(d))
        .attr('x2', d => this.getXFromSeqPosition(d))
        .attr('y1', 0)
        .attr('y2', this._height);
    }
    this._updateHighlight();
  }
}

export default SqrutinyMatrix;
