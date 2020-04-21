var SqrutinyMatrix=function(t){var e={};function n(r){if(e[r])return e[r].exports;var o=e[r]={i:r,l:!1,exports:{}};return t[r].call(o.exports,o,o.exports,n),o.l=!0,o.exports}return n.m=t,n.c=e,n.d=function(t,e,r){n.o(t,e)||Object.defineProperty(t,e,{enumerable:!0,get:r})},n.r=function(t){"undefined"!=typeof Symbol&&Symbol.toStringTag&&Object.defineProperty(t,Symbol.toStringTag,{value:"Module"}),Object.defineProperty(t,"__esModule",{value:!0})},n.t=function(t,e){if(1&e&&(t=n(t)),8&e)return t;if(4&e&&"object"==typeof t&&t&&t.__esModule)return t;var r=Object.create(null);if(n.r(r),Object.defineProperty(r,"default",{enumerable:!0,value:t}),2&e&&"string"!=typeof t)for(var o in t)n.d(r,o,function(e){return t[e]}.bind(null,o));return r},n.n=function(t){var e=t&&t.__esModule?function(){return t.default}:function(){return t};return n.d(e,"a",e),e},n.o=function(t,e){return Object.prototype.hasOwnProperty.call(t,e)},n.p="",n(n.s=8)}([function(t,e){t.exports=d3},function(t,e){function n(e){return t.exports=n=Object.setPrototypeOf?Object.getPrototypeOf:function(t){return t.__proto__||Object.getPrototypeOf(t)},n(e)}t.exports=n},function(t,e){t.exports=function(t,e){if(!(t instanceof e))throw new TypeError("Cannot call a class as a function")}},function(t,e){function n(t,e){for(var n=0;n<e.length;n++){var r=e[n];r.enumerable=r.enumerable||!1,r.configurable=!0,"value"in r&&(r.writable=!0),Object.defineProperty(t,r.key,r)}}t.exports=function(t,e,r){return e&&n(t.prototype,e),r&&n(t,r),t}},function(t,e,n){var r=n(9),o=n(10);t.exports=function(t,e){return!e||"object"!==r(e)&&"function"!=typeof e?o(t):e}},function(t,e,n){var r=n(11);function o(e,n,i){return"undefined"!=typeof Reflect&&Reflect.get?t.exports=o=Reflect.get:t.exports=o=function(t,e,n){var o=r(t,e);if(o){var i=Object.getOwnPropertyDescriptor(o,e);return i.get?i.get.call(n):i.value}},o(e,n,i||e)}t.exports=o},function(t,e,n){var r=n(12);t.exports=function(t,e){if("function"!=typeof e&&null!==e)throw new TypeError("Super expression must either be null or a function");t.prototype=Object.create(e&&e.prototype,{constructor:{value:t,writable:!0,configurable:!0}}),e&&r(t,e)}},function(t,e){t.exports=ProtvistaTrack},function(t,e,n){t.exports=n(13)},function(t,e){function n(t){return(n="function"==typeof Symbol&&"symbol"==typeof Symbol.iterator?function(t){return typeof t}:function(t){return t&&"function"==typeof Symbol&&t.constructor===Symbol&&t!==Symbol.prototype?"symbol":typeof t})(t)}function r(e){return"function"==typeof Symbol&&"symbol"===n(Symbol.iterator)?t.exports=r=function(t){return n(t)}:t.exports=r=function(t){return t&&"function"==typeof Symbol&&t.constructor===Symbol&&t!==Symbol.prototype?"symbol":n(t)},r(e)}t.exports=r},function(t,e){t.exports=function(t){if(void 0===t)throw new ReferenceError("this hasn't been initialised - super() hasn't been called");return t}},function(t,e,n){var r=n(1);t.exports=function(t,e){for(;!Object.prototype.hasOwnProperty.call(t,e)&&null!==(t=r(t)););return t}},function(t,e){function n(e,r){return t.exports=n=Object.setPrototypeOf||function(t,e){return t.__proto__=e,t},n(e,r)}t.exports=n},function(t,e,n){"use strict";n.r(e);var r=n(2),o=n.n(r),i=n(3),s=n.n(i),c=n(4),u=n.n(c),a=n(1),f=n.n(a),l=n(5),h=n.n(l),p=n(6),y=n.n(p),d=n(7),g=n.n(d),b=n(0),x=function(t){function e(){var t;return o()(this,e),(t=u()(this,f()(e).call(this)))._line=Object(b.line)().defined(function(t){return!isNaN(t.y)}).x(function(e){return t.getMiddleXPosition(e.x)}).y(function(e){return t._yScale(e.y)}),t}return y()(e,t),s()(e,[{key:"connectedCallback",value:function(){h()(f()(e.prototype),"connectedCallback",this).call(this),this._height=parseInt(this.getAttribute("height"))||40,this._color=this.getAttribute("color")||"darkgrey",this._cutoffs=this.getAttribute("cutoffs")||[],this._yScale=Object(b.scaleLinear)(),this._xExtent,this._yExtent}},{key:"updateYScale",value:function(){this._xExtent=Object(b.extent)(this._data,function(t){return parseInt(t.x)}),this._yExtent=Object(b.extent)(this._data,function(t){return t.y}),this.xScale.domain(this._xExtent).range([0,this.width]),this._yScale.domain(this._yExtent).range([this._height-this.margin.bottom,this.margin.top])}},{key:"_createTrack",value:function(){this.svg||(this.svg=Object(b.select)(this).append("div").append("svg").attr("width",this.getWidthWithMargins()).attr("height",this._height),this.cutoff_g=this.svg.append("g").attr("class","cutoff"),this.focus=this.svg.append("g").attr("class","focus").style("display","none"),this.focus.append("circle").style("fill","none").attr("stroke","steelblue").attr("r",5),this.focus.append("rect").attr("class","tooltip").style("fill","none").attr("width",50).attr("height",38).attr("x",10).attr("y",-22).attr("rx",4).attr("ry",4),this.focus.append("text").attr("class","value").attr("x",20),this.trackHighlighter.appendHighlightTo(this.svg)),this.updateYScale(),this.refresh()}},{key:"refresh",value:function(){var t=this;this.svg&&(this.svg.selectAll("path").remove(),this.svg.selectAll("line").remove(),this.svg.selectAll("rect").remove(),this.svg.append("path").datum(this._data).attr("d",this._line).attr("fill","none").attr("stroke",this._color).attr("stroke-width","1.2px"),this._cutoffs&&this._cutoffs.length&&this.cutoff_g.selectAll("g").data(this._cutoffs).enter().append("line").style("stroke","rgb(0, 0, 255)").style("stroke-dasharray","5px").style("stroke-width","1.5px").attr("x1",function(e){return t.getMiddleXPosition(e)}).attr("x2",function(e){return t.getMiddleXPosition(e)}).attr("y1",0).attr("y2",this._height),this._updateHighlight())}},{key:"getMiddleXPosition",value:function(t){return(this.getXFromSeqPosition(t)+this.getXFromSeqPosition(t+1))/2}},{key:"data",set:function(t){this._data=void 0,!Array.isArray(t)||t.length<1||(this._data=t.map(function(t){return{x:t.pos,y:t.score}}),this.bisectDate=Object(b.bisector)(function(t){return t.x}).left,this._createTrack())}},{key:"color",set:function(t){this._color=t}},{key:"cutoffs",set:function(t){(!Array.isArray(t)||t.length<1)&&(this._cutoffs=[]),this._cutoffs=t,this.refresh()}}]),e}(g.a);window.customElements&&customElements.define("sqrutiny-matrix",x);e.default=x}]);
//# sourceMappingURL=sqrutiny-matrix.js.map