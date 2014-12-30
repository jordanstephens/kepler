require "kepler/universal_formulation"
require "kepler/lagrange"
require "kepler/laguerre"
require "kepler/param_helper"

module Kepler
  class Orbit
    extend ParamHelper

    attr_accessor :r, :v

    def initialize(r, v)
      @r = r
      @v = v
    end

    def self.from_params(params)
      params = expanded_params(params)

      r_p = perifocal_position(params[:angular_momentum], params[:eccentricity], params[:true_anomaly])
      v_p = perifocal_velocity(params[:angular_momentum], params[:eccentricity], params[:true_anomaly])

      q = transform_matrix(params[:argument_of_periapsis], params[:inclination], params[:right_ascension])

      r = q * r_p
      v = q * v_p

      self.new(r, v)
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
      Math.acos(K.inner_product(h) / h.magnitude).to_deg
    end
    alias_method :i, :inclination

    def node_line
      K.cross_product(h)
    end
    alias_method :n, :node_line

    def right_ascension
      return 0 if n.magnitude == 0
      omega = Math.acos(n[0] / n.magnitude).to_deg
      n[1] < 0 ? (360 - omega) : omega
    end
    alias_method :omega, :right_ascension

    def argument_of_periapsis
      w = Math.acos(n.inner_product(e) / (n.magnitude * e.magnitude)).to_deg
      n[2] < 0 ? (360 - w) : w
    end
    alias_method :w, :argument_of_periapsis

    def true_anomaly
      Math.acos(e.inner_product(r) / (e.magnitude * r.magnitude)).to_deg
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

    # Find the Universal Anomaly by using Laguerre's method to find
    # roots to Kepler's equation in terms of x at time dt
    def universal_anomaly(dt)
      # initial guess of x
      x = Math.sqrt(MU) * (dt / a);

      f = Proc.new { |x| UniversalFormulation.method(:f).call(x, a, @r, @v, dt) }
      df = Proc.new { |x| UniversalFormulation.method(:dfdt).call(x, a, @r, @v) }
      d2f = Proc.new { |x| UniversalFormulation.method(:d2fdt).call(x, a, @r, @v) }

      Laguerre.solve(x, f, df, d2f)
    end

    def update!(dt)
      x = universal_anomaly(dt)
      z = UniversalFormulation.Z(x, a)
      r0 = @r
      v0 = @v

      # make sure you use the same `z` for calculating `@r` and `@v`.
      # this can be tricky because `z` depends on `@r` via `a` so we
      # must be careful to not recalculate `z` between updating `@r`
      # and updating `@v`.
      @r = (Lagrange.f(x, z, @r) * r0) + (Lagrange.g(x, z, dt) * v0)
      @v = (Lagrange.df(x, z, @r, r0) * r0) + (Lagrange.dg(x, z, @r) * v0)

      self
    end

    def update(dt)
      self.class.new(@r, @v).update!(dt)
    end
  end
end
