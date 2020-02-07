var ProtvistaVariationGraph = function(t) {
  var e = {};

  function n(r) {
    if (e[r]) return e[r].exports;
    var o = e[r] = {
      i: r,
      l: !1,
      exports: {}
    };
    return t[r].call(o.exports, o, o.exports, n), o.l = !0, o.exports
  }
  return n.m = t, n.c = e, n.d = function(t, e, r) {
    n.o(t, e) || Object.defineProperty(t, e, {
      enumerable: !0,
      get: r
    })
  }, n.r = function(t) {
    "undefined" != typeof Symbol && Symbol.toStringTag && Object.defineProperty(t, Symbol.toStringTag, {
      value: "Module"
    }), Object.defineProperty(t, "__esModule", {
      value: !0
    })
  }, n.t = function(t, e) {
    if (1 & e && (t = n(t)), 8 & e) return t;
    if (4 & e && "object" == typeof t && t && t.__esModule) return t;
    var r = Object.create(null);
    if (n.r(r), Object.defineProperty(r, "default", {
        enumerable: !0,
        value: t
      }), 2 & e && "string" != typeof t)
      for (var o in t) n.d(r, o, function(e) {
        return t[e]
      }.bind(null, o));
    return r
  }, n.n = function(t) {
    var e = t && t.__esModule ? function() {
      return t.default
    } : function() {
      return t
    };
    return n.d(e, "a", e), e
  }, n.o = function(t, e) {
    return Object.prototype.hasOwnProperty.call(t, e)
  }, n.p = "", n(n.s = 8)
}([function(t, e) {
  t.exports = d3
}, function(t, e) {
  function n(e) {
    return t.exports = n = Object.setPrototypeOf ? Object.getPrototypeOf : function(t) {
      return t.__proto__ || Object.getPrototypeOf(t)
    }, n(e)
  }
  t.exports = n
}, function(t, e) {
  t.exports = function(t, e) {
    if (!(t instanceof e)) throw new TypeError("Cannot call a class as a function")
  }
}, function(t, e) {
  function n(t, e) {
    for (var n = 0; n < e.length; n++) {
      var r = e[n];
      r.enumerable = r.enumerable || !1, r.configurable = !0, "value" in r && (r.writable = !0), Object.defineProperty(t, r.key, r)
    }
  }
  t.exports = function(t, e, r) {
    return e && n(t.prototype, e), r && n(t, r), t
  }
}, function(t, e, n) {
  var r = n(9),
    o = n(10);
  t.exports = function(t, e) {
    return !e || "object" !== r(e) && "function" != typeof e ? o(t) : e
  }
}, function(t, e, n) {
  var r = n(11);

  function o(e, n, i) {
    return "undefined" != typeof Reflect && Reflect.get ? t.exports = o = Reflect.get : t.exports = o = function(t, e, n) {
      var o = r(t, e);
      if (o) {
        var i = Object.getOwnPropertyDescriptor(o, e);
        return i.get ? i.get.call(n) : i.value
      }
    }, o(e, n, i || e)
  }
  t.exports = o
}, function(t, e, n) {
  var r = n(12);
  t.exports = function(t, e) {
    if ("function" != typeof e && null !== e) throw new TypeError("Super expression must either be null or a function");
    t.prototype = Object.create(e && e.prototype, {
      constructor: {
        value: t,
        writable: !0,
        configurable: !0
      }
    }), e && r(t, e)
  }
}, function(t, e) {
  t.exports = ProtvistaTrack
}, function(t, e, n) {
  t.exports = n(13)
}, function(t, e) {
  function n(t) {
    return (n = "function" == typeof Symbol && "symbol" == typeof Symbol.iterator ? function(t) {
      return typeof t
    } : function(t) {
      return t && "function" == typeof Symbol && t.constructor === Symbol && t !== Symbol.prototype ? "symbol" : typeof t
    })(t)
  }

  function r(e) {
    return "function" == typeof Symbol && "symbol" === n(Symbol.iterator) ? t.exports = r = function(t) {
      return n(t)
    } : t.exports = r = function(t) {
      return t && "function" == typeof Symbol && t.constructor === Symbol && t !== Symbol.prototype ? "symbol" : n(t)
    }, r(e)
  }
  t.exports = r
}, function(t, e) {
  t.exports = function(t) {
    if (void 0 === t) throw new ReferenceError("this hasn't been initialised - super() hasn't been called");
    return t
  }
}, function(t, e, n) {
  var r = n(1);
  t.exports = function(t, e) {
    for (; !Object.prototype.hasOwnProperty.call(t, e) && null !== (t = r(t)););
    return t
  }
}, function(t, e) {
  function n(e, r) {
    return t.exports = n = Object.setPrototypeOf || function(t, e) {
      return t.__proto__ = e, t
    }, n(e, r)
  }
  t.exports = n
}, function(t, e, n) {
  "use strict";
  n.r(e);
  var r = n(2),
    o = n.n(r),
    i = n(3),
    a = n.n(i),
    s = n(4),
    c = n.n(s),
    u = n(1),
    f = n.n(u),
    l = n(5),
    p = n.n(l),
    h = n(6),
    y = n.n(h),
    d = n(7),
    b = n.n(d),
    _ = n(0),
    v = function(t) {
      function e() {
        var t;
        return o()(this, e), (t = c()(this, f()(e).call(this)))._line = Object(_.line)().x(function(e) {
          return t.getXFromSeqPosition(e.x)
        }).y(function(e) {
          return t._yScale(e.y)
        }), t
      }
      return y()(e, t), a()(e, [{
        key: "init",
        value: function() {
          this._totals_dataset = {}, this._totals_feature = void 0
        }
      }, {
        key: "connectedCallback",
        value: function() {
          p()(f()(e.prototype), "connectedCallback", this).call(this), this._data = void 0, this._height = parseInt(this.getAttribute("height")) || 40, this._yScale = Object(_.scaleLinear)(), this._xExtent, this._yExtent, this.init()
        }
      }, {
        key: "_createTrack",
        value: function() {
          Object(_.select)(this).selectAll("svg").remove(), this.svg = Object(_.select)(this).append("svg").attr("width", this.width).attr("height", this._height), this.trackHighlighter.appendHighlightTo(this.svg), this._createFeatures(), this.refresh()
        }
      }, {
        key: "_createFeatures",
        value: function() {
          this._xExtent = Object(_.extent)(this._totals_dataset, function(t) {
            return parseInt(t.x)
          }), this._yExtent = Object(_.extent)(this._totals_dataset, function(t) {
            return t.y
          }), this._yExtent[1] += 2, this.xScale.domain(this._xExtent).range([0, this._width]), this._yScale.domain(this._yExtent).range([this._height, 0])
        }
      }, {
        key: "refresh",
        value: function() {
          this.svg && (this.svg.selectAll("path").remove(), this._totals_feature = this.svg.append("path").attr("d", this._line(this._totals_dataset)).attr("fill", "none").attr("stroke", "darkgrey").attr("stroke-width", "1.7px").attr("stroke-dasharray", ".5").attr("transform", "translate(0,0)"), this._updateHighlight())
        }
      }, {
        key: "data",
        set: function(t) {
          if (this._data = t, this.init(), !(this._data.length <= 0)) {
            var e = {},
              n = {};
            this._data.forEach(function(t) {
              void 0 === e[t.start] && (e[t.start] = 0), void 0 === n[t.start] && (n[t.start] = 0), e[t.start] = t.score, void 0 !== t.association && t.association.forEach(function(e) {
                !0 === e.disease && (n[t.start] = t.score)
              })
            }), this._totals_dataset = Object.keys(e).map(function(t) {
              return {
                x: t,
                y: e[t]
              }
            }), this._createTrack()
          }
        }
      }]), e
    }(b.a),
    x = function() {
      customElements.define("protvista-variation-graph", v)
    };
  window.customElements ? x() : document.addEventListener("WebComponentsReady", function() {
    x()
  })
}]);
//# sourceMappingURL=protvista-variation-graph.js.map
