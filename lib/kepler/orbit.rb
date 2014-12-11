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

    def self.from_params(params)
      params = expanded_params(params)

      r_p = perifocal_position(params[:angular_momentum], params[:eccentricity], params[:true_anomaly])
      v_p = perifocal_velocity(params[:angular_momentum], params[:eccentricity], params[:true_anomaly])

      q = transform_matrix(params[:argument_of_periapsis], params[:inclination], params[:right_ascension])

      r = q * r_p
      v = q * v_p

      self.new(r, v)
    end

    def self.expanded_params(params)
      required_params = [
        %i(semimajor_axis eccentricity),
        %i(semilatus_rectum eccentricity),
        %i(apogee perigee)
      ]

      param_set = required_params.find { |rp| rp & params.keys == rp }
      raise ArgumentError, "Invalid parameter set" if param_set.nil?

      default_params = {
        inclination: 0,
        argument_of_periapsis: 0,
        right_ascension: 0,
        true_anomaly: 0,
        eccentricity: 0
      }

      params = default_params.merge(params)

      if param_set.include?(:apogee) && param_set.include?(:perigee)
        params[:semimajor_axis] = ((EARTH_RADIUS * 2) + params[:apogee] + params[:perigee]) / 2
        params[:eccentricity] = (params[:semimajor_axis] / (EARTH_RADIUS + params[:perigee])) - 1
      elsif param_set.include?(:semilatus_rectum)
        params[:semimajor_axis] = params[:semilatus_rectum] / (1 - (params[:eccentricity] ** 2))
      end

      unless param_set.include?(:semilatus_rectum)
        params[:semilatus_rectum] = params[:semimajor_axis] * (1 - (params[:eccentricity] ** 2))
      end

      params[:angular_momentum] = Math.sqrt(params[:semilatus_rectum] * MU)

      params
    end

    # matrix to transform vectors with perifocal basis to vectors with
    # geocentric equatorial basis
    # all args are scalar
    def self.transform_matrix(argument_of_periapsis, inclination, right_ascension)
      w = argument_of_periapsis.to_rad
      i = inclination.to_rad
      omega = right_ascension.to_rad

      Matrix[
        [-Math.sin(omega) * Math.cos(i) * Math.sin(w) + (Math.cos(omega) * Math.cos(w)),
         -Math.sin(omega) * Math.cos(i) * Math.cos(w) - (Math.cos(omega) * Math.sin(w)),
         Math.sin(omega) * Math.sin(i)],
        [Math.cos(omega) * Math.cos(i) * Math.sin(w) + (Math.sin(omega) * Math.cos(w)),
         Math.cos(omega) * Math.cos(i) * Math.cos(w) - (Math.sin(omega) * Math.sin(w)),
         -Math.cos(omega) * Math.sin(i)],
        [Math.sin(i) * Math.sin(w),
         Math.sin(i) * Math.cos(w),
         Math.cos(i)]
      ]
    end

    # all args are scalar
    def self.perifocal_position(angular_momentum, eccentricity, true_anomaly)
      h = angular_momentum
      e = eccentricity
      theta = true_anomaly

      ((h ** 2) / MU) *
      (1 / (1 + (e * Math.cos(theta.to_rad)))) *
      Vector[Math.cos(theta.to_rad), Math.sin(theta.to_rad), 0]
    end

    # all args are scalar
    def self.perifocal_velocity(angular_momentum, eccentricity, true_anomaly)
      h = angular_momentum
      e = eccentricity
      theta = true_anomaly

      (MU / h) *
      Vector[-Math.sin(theta.to_rad), e + Math.cos(theta.to_rad), 0]
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
  end
end
