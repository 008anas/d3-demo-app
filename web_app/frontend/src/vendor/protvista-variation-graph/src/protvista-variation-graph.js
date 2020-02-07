import ProtvistaTrack from "protvista-track";
import {
  scaleLinear,
  select,
  event as d3Event,
  line,
  extent,
  curveBasis
} from "d3";

class ProtvistaVariationGraph extends ProtvistaTrack {
  constructor() {
    super();
    this._line = line()
      .x(d => this.getXFromSeqPosition(d.x))
      .y(d => this._yScale(d.y));
  }

  init() {
    this._totals_dataset = {};
    this._totals_feature = undefined;
  }

  connectedCallback() {
    super.connectedCallback();

    this._data = undefined;

    this._height = parseInt(this.getAttribute("height")) || 40;
    this._yScale = scaleLinear();
    this._xExtent;
    this._yExtent;
    this.init();
  }

  set data(data) {
    this._data = data;
    this.init();
    if (this._data.variants.length <= 0) {
      return;
    }

    let totalMap = {};
    let diseaseMap = {};

    this._data.variants.forEach(v => {
      if ("undefined" === typeof totalMap[v.start]) {
        totalMap[v.start] = 0;
      }
      if ("undefined" === typeof diseaseMap[v.start]) {
        diseaseMap[v.start] = 0;
      }
      totalMap[v.start] = v.score;
      if ("undefined" !== typeof v.association) {
        v.association.forEach(a => {
          if (true === a.disease) {
            diseaseMap[v.start] = v.score;
          }
        });
      }
    });
    this._totals_dataset = Object.keys(totalMap).map(d => {
      return {
        x: d,
        y: totalMap[d]
      };
    });
    this._createTrack();
  }

  _createTrack() {
    select(this)
      .selectAll("svg")
      .remove();
    this.svg = select(this)
      .append("svg")
      .attr("width", this.width)
      .attr("height", this._height);
    this.trackHighlighter.appendHighlightTo(this.svg);
    // Create the visualisation here
    this._createFeatures();
    this.refresh();
  }

  _createFeatures() {
    this._xExtent = extent(this._totals_dataset, d => parseInt(d.x));
    this._yExtent = extent(this._totals_dataset, d => d.y);

    // just a bit of padding on the top
    this._yExtent[1] += 2;

    this.xScale.domain(this._xExtent).range([0, this._width]);
    this._yScale.domain(this._yExtent).range([this._height, 0]);
  }

  refresh() {
    if (!this.svg) return;
    this.svg.selectAll("path").remove();
    this._totals_feature = this.svg
      .append("path")
      .attr("d", this._line(this._totals_dataset))
      .attr("fill", "none")
      .attr("stroke", "darkgrey")
      .attr("stroke-width", "1px")
      .attr("stroke-dasharray", ".5")
      .attr("transform", "translate(0,0)");
    this._updateHighlight();
  }
}

export default ProtvistaVariationGraph;
