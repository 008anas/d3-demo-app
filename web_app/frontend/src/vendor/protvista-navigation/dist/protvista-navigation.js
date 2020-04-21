var ProtvistaNavigation=function(t){var e={};function i(n){if(e[n])return e[n].exports;var r=e[n]={i:n,l:!1,exports:{}};return t[n].call(r.exports,r,r.exports,i),r.l=!0,r.exports}return i.m=t,i.c=e,i.d=function(t,e,n){i.o(t,e)||Object.defineProperty(t,e,{enumerable:!0,get:n})},i.r=function(t){"undefined"!=typeof Symbol&&Symbol.toStringTag&&Object.defineProperty(t,Symbol.toStringTag,{value:"Module"}),Object.defineProperty(t,"__esModule",{value:!0})},i.t=function(t,e){if(1&e&&(t=i(t)),8&e)return t;if(4&e&&"object"==typeof t&&t&&t.__esModule)return t;var n=Object.create(null);if(i.r(n),Object.defineProperty(n,"default",{enumerable:!0,value:t}),2&e&&"string"!=typeof t)for(var r in t)i.d(n,r,function(e){return t[e]}.bind(null,r));return n},i.n=function(t){var e=t&&t.__esModule?function(){return t.default}:function(){return t};return i.d(e,"a",e),e},i.o=function(t,e){return Object.prototype.hasOwnProperty.call(t,e)},i.p="",i(i.s=8)}([function(t,e){t.exports=d3},function(t,e){function i(e,n){return t.exports=i=Object.setPrototypeOf||function(t,e){return t.__proto__=e,t},i(e,n)}t.exports=i},function(t,e){function i(e){return t.exports=i=Object.setPrototypeOf?Object.getPrototypeOf:function(t){return t.__proto__||Object.getPrototypeOf(t)},i(e)}t.exports=i},function(t,e){t.exports=function(t,e){if(!(t instanceof e))throw new TypeError("Cannot call a class as a function")}},function(t,e){function i(t,e){for(var i=0;i<e.length;i++){var n=e[i];n.enumerable=n.enumerable||!1,n.configurable=!0,"value"in n&&(n.writable=!0),Object.defineProperty(t,n.key,n)}}t.exports=function(t,e,n){return e&&i(t.prototype,e),n&&i(t,n),t}},function(t,e,i){var n=i(9),r=i(10);t.exports=function(t,e){return!e||"object"!==n(e)&&"function"!=typeof e?r(t):e}},function(t,e,i){var n=i(1);t.exports=function(t,e){if("function"!=typeof e&&null!==e)throw new TypeError("Super expression must either be null or a function");t.prototype=Object.create(e&&e.prototype,{constructor:{value:t,writable:!0,configurable:!0}}),e&&n(t,e)}},function(t,e,i){var n=i(2),r=i(1),s=i(11),o=i(12);function a(e){var i="function"==typeof Map?new Map:void 0;return t.exports=a=function(t){if(null===t||!s(t))return t;if("function"!=typeof t)throw new TypeError("Super expression must either be null or a function");if(void 0!==i){if(i.has(t))return i.get(t);i.set(t,e)}function e(){return o(t,arguments,n(this).constructor)}return e.prototype=Object.create(t.prototype,{constructor:{value:e,enumerable:!1,writable:!0,configurable:!0}}),r(e,t)},a(e)}t.exports=a},function(t,e,i){t.exports=i(13)},function(t,e){function i(t){return(i="function"==typeof Symbol&&"symbol"==typeof Symbol.iterator?function(t){return typeof t}:function(t){return t&&"function"==typeof Symbol&&t.constructor===Symbol&&t!==Symbol.prototype?"symbol":typeof t})(t)}function n(e){return"function"==typeof Symbol&&"symbol"===i(Symbol.iterator)?t.exports=n=function(t){return i(t)}:t.exports=n=function(t){return t&&"function"==typeof Symbol&&t.constructor===Symbol&&t!==Symbol.prototype?"symbol":i(t)},n(e)}t.exports=n},function(t,e){t.exports=function(t){if(void 0===t)throw new ReferenceError("this hasn't been initialised - super() hasn't been called");return t}},function(t,e){t.exports=function(t){return-1!==Function.toString.call(t).indexOf("[native code]")}},function(t,e,i){var n=i(1);function r(e,i,s){return!function(){if("undefined"==typeof Reflect||!Reflect.construct)return!1;if(Reflect.construct.sham)return!1;if("function"==typeof Proxy)return!0;try{return Date.prototype.toString.call(Reflect.construct(Date,[],function(){})),!0}catch(t){return!1}}()?t.exports=r=function(t,e,i){var r=[null];r.push.apply(r,e);var s=new(Function.bind.apply(t,r));return i&&n(s,i.prototype),s}:t.exports=r=Reflect.construct,r.apply(null,arguments)}t.exports=r},function(t,e,i){"use strict";i.r(e);var n=i(3),r=i.n(n),s=i(4),o=i.n(s),a=i(5),h=i.n(a),u=i(2),l=i.n(u),c=i(6),p=i.n(c),d=i(7),f=i.n(d),_=i(0),y=function(t){function e(){var t;return r()(this,e),(t=h()(this,l()(e).call(this)))._x=null,t._padding=0,t.dontDispatch=!1,t}return p()(e,t),o()(e,[{key:"_refreshWidth",value:function(){this.style.display="block",this.style.width="100%",this.width=this.offsetWidth,this.width>0&&(this._padding=10)}},{key:"connectedCallback",value:function(){this._refreshWidth(),this.closest("protvista-manager")&&(this.manager=this.closest("protvista-manager"),this.manager.register(this)),this._displaystart=parseFloat(this.getAttribute("displaystart"))||1,this._highlightStart=parseFloat(this.getAttribute("highlightStart")),this._highlightEnd=parseFloat(this.getAttribute("highlightEnd")),this._onResize=this._onResize.bind(this),this._length=parseFloat(this.getAttribute("length")),this.init()}},{key:"init",value:function(){this._displayend=parseFloat(this.getAttribute("displayend"))||this._length,this._length>0&&this._createNavRuler()}},{key:"disconnectedCallback",value:function(){this.manager&&this.manager.unregister(this),this._ro&&this._ro.unobserve(this),window.removeEventListener("resize",this._onResize)}},{key:"attributeChangedCallback",value:function(t,e,i){e!==i&&(this["_".concat(t)]=parseFloat(i),this._updateNavRuler())}},{key:"_createNavRuler",value:function(){var t=this;this._x=Object(_.scaleLinear)().range([this._padding,this.width-this._padding]),this._x.domain([1,this._length]),this._svg=Object(_.select)(this).append("div").attr("class","navigator").append("svg").attr("id","nav-svg").attr("width",this.width).attr("height",40),this._xAxis=Object(_.axisBottom)(this._x),this._displaystartLabel=this._svg.append("text").attr("class","start-label").attr("x",0).attr("y",40-this._padding),this._displayendLabel=this._svg.append("text").attr("class","end-label").attr("x",this.width).attr("y",40-this._padding).attr("text-anchor","end"),this._axis=this._svg.append("g").attr("class","x axis").call(this._xAxis),this._viewport=Object(_.brushX)().extent([[this._padding,0],[this.width-this._padding,20.4]]).on("brush",function(){_.event.selection&&(t._displaystart=Object(_.format)("d")(t._x.invert(_.event.selection[0])),t._displayend=Object(_.format)("d")(t._x.invert(_.event.selection[1])),t.dontDispatch||t.dispatchEvent(new CustomEvent("change",{detail:{displayend:t._displayend,displaystart:t._displaystart,extra:{transform:_.event.transform}},bubbles:!0,cancelable:!0})),t._updateLabels(),t._updatePolygon())}),this._brushG=this._svg.append("g").attr("class","brush").call(this._viewport),this._brushG.call(this._viewport.move,[this._x(this._displaystart),this._x(this._displayend)]),this.polygon=this._svg.append("polygon").attr("class","zoom-polygon").attr("fill","#777").attr("fill-opacity","0.3"),this._updateNavRuler(),"ResizeObserver"in window&&(this._ro=new ResizeObserver(this._onResize),this._ro.observe(this)),window.addEventListener("resize",this._onResize)}},{key:"_onResize",value:function(){this._refreshWidth(),this._x=this._x.range([this._padding,this.width-this._padding]),this._svg.attr("width",this.width),this._viewport.extent([[this._padding,0],[this.width-this._padding,20.4]]),this._brushG.call(this._viewport),this._updateNavRuler()}},{key:"_updateNavRuler",value:function(){this._x&&(this._x.domain([1,this._length]),this._axis.call(this._xAxis),this._updatePolygon(),this._updateLabels(),this._brushG&&(this.dontDispatch=!0,this._brushG.call(this._viewport.move,[this._x(this._displaystart),this._x(this._displayend)]),this.dontDispatch=!1))}},{key:"_updateLabels",value:function(){this._displaystartLabel&&this._displaystartLabel.text(this._displaystart),this._displayendLabel&&this._displayendLabel.attr("x",this.width).text(this._displayend)}},{key:"_updatePolygon",value:function(){this.polygon&&this.polygon.attr("points","".concat(this._x(this._displaystart),",").concat(20,"\n        ").concat(this._x(this._displayend),",").concat(20,"\n        ").concat(this.width,",").concat(40,"\n        0,").concat(40))}},{key:"width",get:function(){return this._width},set:function(t){this._width=t}},{key:"length",set:function(t){this._length=t,this.init()}}],[{key:"observedAttributes",get:function(){return["length","displaystart","displayend","highlightStart","highlightEnd","width"]}}]),e}(f()(HTMLElement)),b=function(){customElements.define("protvista-navigation",y)};window.customElements?b():document.addEventListener("WebComponentsReady",function(){b()})}]);
//# sourceMappingURL=protvista-navigation.js.map