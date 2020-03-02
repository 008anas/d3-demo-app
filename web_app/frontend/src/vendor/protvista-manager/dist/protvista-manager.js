var ProtvistaManager=function(t){var e={};function n(r){if(e[r])return e[r].exports;var o=e[r]={i:r,l:!1,exports:{}};return t[r].call(o.exports,o,o.exports,n),o.l=!0,o.exports}return n.m=t,n.c=e,n.d=function(t,e,r){n.o(t,e)||Object.defineProperty(t,e,{enumerable:!0,get:r})},n.r=function(t){"undefined"!=typeof Symbol&&Symbol.toStringTag&&Object.defineProperty(t,Symbol.toStringTag,{value:"Module"}),Object.defineProperty(t,"__esModule",{value:!0})},n.t=function(t,e){if(1&e&&(t=n(t)),8&e)return t;if(4&e&&"object"==typeof t&&t&&t.__esModule)return t;var r=Object.create(null);if(n.r(r),Object.defineProperty(r,"default",{enumerable:!0,value:t}),2&e&&"string"!=typeof t)for(var o in t)n.d(r,o,function(e){return t[e]}.bind(null,o));return r},n.n=function(t){var e=t&&t.__esModule?function(){return t.default}:function(){return t};return n.d(e,"a",e),e},n.o=function(t,e){return Object.prototype.hasOwnProperty.call(t,e)},n.p="",n(n.s=7)}([function(t,e){function n(e,r){return t.exports=n=Object.setPrototypeOf||function(t,e){return t.__proto__=e,t},n(e,r)}t.exports=n},function(t,e){function n(e){return t.exports=n=Object.setPrototypeOf?Object.getPrototypeOf:function(t){return t.__proto__||Object.getPrototypeOf(t)},n(e)}t.exports=n},function(t,e){t.exports=function(t,e){if(!(t instanceof e))throw new TypeError("Cannot call a class as a function")}},function(t,e){function n(t,e){for(var n=0;n<e.length;n++){var r=e[n];r.enumerable=r.enumerable||!1,r.configurable=!0,"value"in r&&(r.writable=!0),Object.defineProperty(t,r.key,r)}}t.exports=function(t,e,r){return e&&n(t.prototype,e),r&&n(t,r),t}},function(t,e,n){var r=n(8),o=n(9);t.exports=function(t,e){return!e||"object"!==r(e)&&"function"!=typeof e?o(t):e}},function(t,e,n){var r=n(0);t.exports=function(t,e){if("function"!=typeof e&&null!==e)throw new TypeError("Super expression must either be null or a function");t.prototype=Object.create(e&&e.prototype,{constructor:{value:t,writable:!0,configurable:!0}}),e&&r(t,e)}},function(t,e,n){var r=n(1),o=n(0),i=n(10),u=n(11);function a(e){var n="function"==typeof Map?new Map:void 0;return t.exports=a=function(t){if(null===t||!i(t))return t;if("function"!=typeof t)throw new TypeError("Super expression must either be null or a function");if(void 0!==n){if(n.has(t))return n.get(t);n.set(t,e)}function e(){return u(t,arguments,r(this).constructor)}return e.prototype=Object.create(t.prototype,{constructor:{value:e,enumerable:!1,writable:!0,configurable:!0}}),o(e,t)},a(e)}t.exports=a},function(t,e,n){t.exports=n(12)},function(t,e){function n(t){return(n="function"==typeof Symbol&&"symbol"==typeof Symbol.iterator?function(t){return typeof t}:function(t){return t&&"function"==typeof Symbol&&t.constructor===Symbol&&t!==Symbol.prototype?"symbol":typeof t})(t)}function r(e){return"function"==typeof Symbol&&"symbol"===n(Symbol.iterator)?t.exports=r=function(t){return n(t)}:t.exports=r=function(t){return t&&"function"==typeof Symbol&&t.constructor===Symbol&&t!==Symbol.prototype?"symbol":n(t)},r(e)}t.exports=r},function(t,e){t.exports=function(t){if(void 0===t)throw new ReferenceError("this hasn't been initialised - super() hasn't been called");return t}},function(t,e){t.exports=function(t){return-1!==Function.toString.call(t).indexOf("[native code]")}},function(t,e,n){var r=n(0);function o(e,n,i){return!function(){if("undefined"==typeof Reflect||!Reflect.construct)return!1;if(Reflect.construct.sham)return!1;if("function"==typeof Proxy)return!0;try{return Date.prototype.toString.call(Reflect.construct(Date,[],function(){})),!0}catch(t){return!1}}()?t.exports=o=function(t,e,n){var o=[null];o.push.apply(o,e);var i=new(Function.bind.apply(t,o));return n&&r(i,n.prototype),i}:t.exports=o=Reflect.construct,o.apply(null,arguments)}t.exports=o},function(t,e,n){"use strict";n.r(e);var r=n(2),o=n.n(r),i=n(3),u=n.n(i),a=n(4),c=n.n(a),s=n(1),l=n.n(s),f=n(5),p=n.n(f),y=n(6),b=function(t){function e(){var t;return o()(this,e),(t=c()(this,l()(e).call(this))).protvistaElements=new Set,t.attributeValues=new Map,t.propertyValues=new Map,t}return p()(e,t),u()(e,[{key:"attributeChangedCallback",value:function(t,e,n){if(e!==n&&"attributes"===t){if(this._attributes=n.split(" "),-1!==this._attributes.indexOf("type"))throw new Error("'type' can't be used as a protvista attribute");if(-1!==this._attributes.indexOf("value"))throw new Error("'value' can't be used as a protvista attribute")}}},{key:"register",value:function(t){this.protvistaElements.add(t)}},{key:"unregister",value:function(t){this.protvistaElements.delete(t)}},{key:"_polyfillElementClosest",value:function(){var t=this;Element.prototype.matches||(Element.prototype.matches=Element.prototype.msMatchesSelector||Element.prototype.webkitMatchesSelector),Element.prototype.closest||(Element.prototype.closest=function(e){var n=t;do{if(n.matches(e))return n;n=n.parentElement||n.parentNode}while(null!==n&&1===n.nodeType);return null})}},{key:"applyAttributes",value:function(){var t=this;this.protvistaElements.forEach(function(e){t.attributeValues.forEach(function(t,n){!1===t?e.removeAttribute(n):e.setAttribute(n,"boolean"==typeof t?"":t)})})}},{key:"applyProperties",value:function(t){var e=this;if(t){var n=this.querySelector("#".concat(t));this.propertyValues.forEach(function(t,e){n[e]=t})}else this.protvistaElements.forEach(function(t){e.propertyValues.forEach(function(e,n){t[n]=e})})}},{key:"_changeListener",value:function(t){var e=this;if(t&&t.detail)switch(t.detail.handler){case"property":this.propertyValues.set(t.detail.type,t.detail.value),this.applyProperties(t.detail.for);break;default:-1!==this._attributes.indexOf(t.detail.type)&&this.attributeValues.set(t.detail.type,t.detail.value),Object.keys(t.detail).forEach(function(n){-1!==e._attributes.indexOf(n)&&e.attributeValues.set(n,t.detail[n])}),this.applyAttributes()}}},{key:"connectedCallback",value:function(){this.addEventListener("change",this._changeListener)}}],[{key:"observedAttributes",get:function(){return["attributes"]}}]),e}(n.n(y)()(HTMLElement)),d=function(){customElements.define("protvista-manager",b)};window.customElements?d():document.addEventListener("WebComponentsReady",function(){d()})}]);
//# sourceMappingURL=protvista-manager.js.map