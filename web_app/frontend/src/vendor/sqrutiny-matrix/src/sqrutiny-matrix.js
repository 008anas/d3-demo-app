import {
  scaleLinear,
  line,
  curveMonotoneX,
  select,
  max,
  bisector,
  extent,
  mouse,
  format
} from 'd3';

import ProtvistaTrack from 'protvista-track';

class SqrutinyMatrix extends ProtvistaTrack {

  constructor() {
    super();
    this._line = line()
      .defined(d => !isNaN(d.y))
      .x(d => this.getMiddleXPosition(d.x))
      .y(d => this._yScale(d.y));
    // .curve(curveMonotoneX);
  }

  connectedCallback() {
    super.connectedCallback();

    this._height = parseInt(this.getAttribute('height')) || 40;
    this._color = this.getAttribute('color') || 'darkgrey';
    this._cutoffs = this.getAttribute('cutoffs') || [];
    this._yScale = scaleLinear();
    this._xExtent;
    this._yExtent;
  }

  set data(data) {
    this._data = undefined;
    if (!Array.isArray(data) || data.length < 1) {
      return;
    }
    this._data = data.map(d => ({
      x: d.pos,
      y: d.score
    }));
    this.bisectDate = bisector(d => d.x).left;
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

      this.focus = this.svg.append('g')
        .attr('class', 'focus')
        .style('display', 'none');

      this.focus.append('circle')
        .style("fill", "none")
        .attr('stroke', 'steelblue')
        .attr('r', 5);

      this.focus.append('rect')
        .attr('class', 'tooltip')
        .style('fill', 'none')
        .attr('width', 50)
        .attr('height', 38)
        .attr('x', 10)
        .attr('y', -22)
        .attr('rx', 4)
        .attr('ry', 4);

      this.focus.append('text')
        .attr('class', 'value')
        .attr('x', 20);

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
    this.svg.selectAll('rect').remove();
    this.svg.append('path')
      .datum(this._data)
      .attr('d', this._line)
      .attr('fill', 'none')
      .attr('stroke', this._color)
      .attr('stroke-width', `1.2px`);
    // this.svg.append('rect')
    //   .attr('class', 'overlay')
    //   .style('fill', 'none')
    //   .attr('pointer-events', 'all')
    //   .attr('width', this.width - this.margin.right)
    //   .attr('height', this.height)
    //   .attr('transform', `translate(${this.margin.left},0)`)
    //   .on('mouseover', () => this.focus.style('display', null))
    //   .on('mouseout', () => this.focus.style('display', 'none'))
    //   .on('touchmove mousemove', () => {
    //     var x0 = this.xScale.invert(mouse(this)[0] - this.margin.left),
    //       i = this.bisectDate(this._data, x0),
    //       d = this._data[i];
    //     if (!d) {
    //       return;
    //     }
    //     this.focus.attr('transform', `translate(${this.getMiddleXPosition(d.x)},${this._yScale(d.y)})`);
    //     this.focus.select('.value')
    //       .style('font-weight', 'bold')
    //       .text(!(d.y % 1) ? format(',')(d.y) : format(',.2f')(d.y));
    //   });
    if (this._cutoffs && this._cutoffs.length) {
      this.cutoff_g
        .selectAll('g')
        .data(this._cutoffs)
        .enter()
        .append('line')
        .style('stroke', 'rgb(0, 0, 255)')
        .style('stroke-dasharray', '5px')
        .style('stroke-width', '1.5px')
        .attr('x1', d => this.getMiddleXPosition(d))
        .attr('x2', d => this.getMiddleXPosition(d))
        .attr('y1', 0)
        .attr('y2', this._height);
    }
    this._updateHighlight();
  }

  getMiddleXPosition(pos) {
    return (this.getXFromSeqPosition(pos) + this.getXFromSeqPosition(pos + 1)) / 2;
  }
}

export default SqrutinyMatrix;
