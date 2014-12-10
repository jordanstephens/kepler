require "kepler/universal_formulation"
require "kepler/lagrange"
require "kepler/laguerre"

module Kepler
  class Orbit
    include UniversalFormulation
    include Lagrange

    attr_accessor :r, :v

    def initialize(r, v)
      @r = r
      @v = v
    end

    def angular_momentum
      @r.cross_product(@v)
    end
    alias_method :h, :angular_momentum

    def radial_velocity
      @r.inner_product(@v) / @r.magnitude
    end

    def eccentricity
      (1 / MU) * (
        (((@v.magnitude ** 2) - (MU / @r.magnitude)) * @r) -
        (@r.magnitude * radial_velocity) * @v
      )
    end
    alias_method :e, :eccentricity

    def semimajor_axis
      ((h.magnitude ** 2) / MU) * (1 / (1 - (e.magnitude ** 2)))
    end
    alias_method :a, :semimajor_axis

    def semilatus_rectum
      (h.magnitude ** 2) / MU
    end
    alias_method :p, :semilatus_rectum

    def inclination
      to_deg(Math.acos(K.inner_product(h) / h.magnitude))
    end
    alias_method :i, :inclination

    def node_line
      K.cross_product(h)
    end
    alias_method :n, :node_line

    def right_ascension
      return 0 if n.magnitude == 0
      omega = to_deg(Math.acos(n[0] / n.magnitude))
      n[1] < 0 ? (360 - omega) : omega
    end
    alias_method :omega, :right_ascension

    def argument_of_periapsis
      w = to_deg(Math.acos(n.inner_product(e) / (n.magnitude * e.magnitude)))
      n[2] < 0 ? (360 - w) : w
    end
    alias_method :w, :argument_of_periapsis

    def true_anomaly
      to_deg(Math.acos(e.inner_product(r) / (e.magnitude * r.magnitude)))
    end
    alias_method :theta, :true_anomaly

    def periapsis
      ((h.magnitude ** 2) / MU) * (1 / (1 + e.magnitude * Math.cos(0)))
    end
    alias_method :perigee, :periapsis

    def apoapsis
      ((h.magnitude ** 2) / MU) * (1 / (1 + e.magnitude * Math.cos(Math::PI)))
    end
    alias_method :apogee, :apoapsis

    def period
      (2 * Math::PI / Math.sqrt(MU)) * Math.sqrt(a ** 3)
    end
    alias_method :T, :period

    def universal_anomaly(dt)
      # initial guess of x
      x = Math.sqrt(MU) * (dt / a);

      f = Proc.new { |x| method(:uf_F).call(x, dt) }
      df = method(:uf_dFdt)
      d2f = method(:uf_d2Fdt)

      Laguerre.solve(x, f, df, d2f)
    end

    def update!(dt)
      x = universal_anomaly(dt)
      z = Z(x)
      c = C(z)
      s = S(z)
      r0 = @r
      v0 = @v

      # make sure you use the same `z` for calculating `@r` and `@v`.
      # this can be tricky because `z` depends on `@r` via `a` so we
      # must be careful to not recalculate `z` between updating `@r`
      # and updating `@v`.
      @r = (f(x, z) * r0) + (g(x, z, dt) * v0)
      @v = (df(x, z, r0) * r0) + (dg(x, z) * v0)

      self
    end

    def update(dt)
      self.class.new(@r, @v).update!(dt)
    end

    private

    def to_deg(rad)
      rad * 180 / Math::PI
    end
  end
end
