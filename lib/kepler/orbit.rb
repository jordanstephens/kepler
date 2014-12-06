module Kepler
  class Orbit
    include Stumpff

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

      # find root of f with Laguerre's method
      5.times do |i|
        z = Z(x)
        s = S(z)
        c = C(z)

        f = ((1 - (r.magnitude / a)) * s * (x ** 3)) +
            ((r.inner_product(v) / Math.sqrt(MU)) * c * (x ** 2)) +
            (r.magnitude * x) -
            (Math.sqrt(MU) * dt)

        df = (c * (x ** 2)) +
             ((r.inner_product(v) / Math.sqrt(MU)) * (1 - (s * z)) * x) +
             (r.magnitude * (1 - (c * z)))

        ddf = ((1 - (r.magnitude / a)) * (1 - (s * z)) * x) +
              ((r.inner_product(v) / Math.sqrt(MU)) * (1 - (c * z)))

        delta = 2 * Math.sqrt((4 * (df ** 2)) - (5 * f * ddf))
        dx = (5 * f) / (df + ((df.abs / df) * delta))
        x = x - dx
      end

      x
    end
    alias_method :x, :universal_anomaly

    def update!(dt)
      x = universal_anomaly(dt)
      z = Z(x)
      c = C(z)
      s = S(z)
      r0 = @r
      v0 = @v

      f = 1 - ((x ** 2) / r0.magnitude) * c
      g = dt - ((1 / Math.sqrt(MU)) * (x ** 3) * s)

      @r = (f * r0) + (g * v0)

      df = (Math.sqrt(MU) / (@r.magnitude * r0.magnitude)) * (s * z - 1) * x
      dg = (1 - (((x ** 2) / r.magnitude) * c))

      @v = (df * r0) + (dg * v0)

      self
    end

    def update(dt)
      self.class.new(@r, @v).update!(dt)
    end

    def Z(x)
      (x ** 2) / a
    end

    private
    def to_deg(rad)
      rad * 180 / Math::PI
    end
  end
end
