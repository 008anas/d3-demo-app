import {
  scaleLinear,
  line,
  select,
  min,
  max,
  bisector,
  mouse,
  format
} from 'd3';

import ProtvistaZoomable from 'protvista-zoomable';

class SqrutinyMatrix extends ProtvistaZoomable {

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
    this.refresh();
  }

  set cutoffs(cutoffs) {
    if (!Array.isArray(cutoffs) || cutoffs.length < 1) {
      this._cutoffs = [];
    }
    this._cutoffs = cutoffs;
    this.updateCutoffs();
  }

  // eslint-disable-next-line class-methods-use-this
  get margin() {
    return {
      top: 20,
      right: 15,
      bottom: 20,
      left: 10
    };
  }

  constructor() {
    super();
    this._line = line()
      .defined(d => !isNaN(d.y))
      .x(d => this.getMiddleXPos(d.x))
      .y(d => this.yScale(d.y));
  }

  connectedCallback() {
    super.connectedCallback();

    this._color = this.getAttribute('color') || 'darkgrey';
    this._cutoffs = this.getAttribute('cutoffs') || [];
    this.yScale = scaleLinear();
  }

  _createTrack() {
    if (!this.svg) {

      this.svg = select(this)
        .append('div')
        .append('svg')
        .attr('width', this.width)
        .attr('height', this.height);

      this.inner_width = this.svg.attr('width') - this.margin.left;
      this.inner_height = this.svg.attr('height') - this.margin.top - this.margin.bottom;

      this.line_chart = this.svg.append('g')
        .attr('class', 'focus')
        .attr('transform', `translate(${this.margin.left}, ${this.margin.top})`)
        .attr('clip-path', 'url(#clip)');

      this.tooltipT = this.svg.append('g')
        .attr('class', 'tooltip 1')
        .style('display', 'none');

      this.tooltipB = this.svg.append('g')
        .attr('class', 'tooltip 2')
        .style('display', 'none');

      this.tooltipT.append('circle')
        .attr('stroke', 'steelblue')
        .attr('fill', 'none')
        .attr('r', 4);

      this.tooltipT.append('text')
        .attr('class', 'value')
        .style('font', '11px sans-serif')
        .attr('x', -3)
        .attr('y', -7);

      this.tooltipB.append('text')
        .attr('class', 'value')
        .style('font', '10px sans-serif')
        .style('fill', '#ababab')
        .attr('x', -5)
        .attr('y', 15);

      this.cutoff_g = this.svg.append('g')
        .attr('transform', `translate(${this.margin.left}, ${this.margin.top})`);

      this.rect = this.svg.append('rect')
        .attr('class', 'overlay')
        .attr('fill', 'none')
        .attr('width', this.inner_width)
        .attr('height', this.inner_height + this.margin.top)
        .attr('transform', `translate(${this.margin.left},0)`)
        .attr('pointer-events', 'all')
        .on('mouseover', () => this.svg.selectAll('.tooltip').style('display', null))
        .on('mouseout', () => this.svg.selectAll('.tooltip').style('display', 'none'))
        .on('mousemove', () => {
          this.svg.selectAll('.tooltip').style('display', null);
          var x0 = this.xScale.invert(mouse(this.rect.node())[0] - this.getSingleBaseWidth() / 2),
            i = this.bisectDate(this._data, x0),
            d = this._data[i];
          if (!d) {
            return;
          }
          this.tooltipT.attr('transform', `translate(${this.getMiddleXPosMargin(d.x)},${this.yScale(d.y) + this.margin.top})`);
          this.tooltipT.select('.value').text(!(d.y % 1) ? format(',')(d.y) : format(',.2f')(d.y));
          this.tooltipB.attr('transform', `translate(${this.getMiddleXPosMargin(d.x)},${this.inner_height + this.margin.top})`);
          this.tooltipB.select('.value').text(d.x);
        });

      this.trackHighlighter.appendHighlightTo(this.svg);
    }

    // Create the visualisation here
    this.updateYScale();
    this.refresh();
  }

  updateYScale() {
    this.yScale.domain([min(this._data.map(d => d.y)), max(this._data.map(d => d.y))])
      .range([this.inner_height, 0]);
  }

  refresh() {
    if (!this.svg) return;

    this.svg.selectAll('path.line').remove();
    this.svg.selectAll('g.tooltip').style('display', 'none');

    this.line_chart.append('path')
      .datum(this._data)
      .attr('class', 'line')
      .attr('fill', 'none')
      .attr('stroke', this._color)
      .attr('stroke-width', 1)
      .attr('d', this._line);

    this.updateCutoffs();
  }

  updateCutoffs() {
    if (!this.svg) return;

    this.svg.selectAll('line.cutoff').remove();

    if (this._cutoffs && this._cutoffs.length) {
      this.cutoff_g
        .selectAll('g')
        .data(this._cutoffs)
        .enter()
        .append('line')
        .attr('class', 'cutoff')
        .style('stroke', 'rgb(0, 0, 255)')
        .style('stroke-dasharray', '4px')
        .style('stroke-width', '1.5px')
        .attr('x1', d => this.getMiddleXPos(d))
        .attr('x2', d => this.getMiddleXPos(d))
        .attr('y1', 0)
        .attr('y2', this.inner_height);
    }

    this.trackHighlighter.updateHighlight();
  }

  getMiddleXPosMargin(pos) {
    return this.getXFromSeqPosition(pos) + (this.getSingleBaseWidth() / 2);
  }

  getMiddleXPos(pos) {
    return this.xScale(pos) + (this.getSingleBaseWidth() / 2);
  }
}

export default SqrutinyMatrix;
