
module Kepler
  # Lagrange Coefficients in terms of the change of universal anomaly
  module Lagrange
    class << self
      include Stumpff

      def f(x, z, r)
        1 - ((x ** 2) / r.magnitude) * c(z)
      end

      def g(x, z, mu, dt)
        dt - ((1 / Math.sqrt(mu)) * (x ** 3) * s(z))
      end

      def df(x, z, r, r0, mu)
        (Math.sqrt(mu) / (r.magnitude * r0.magnitude)) * (s(z) * z - 1) * x
      end

      def dg(x, z, r)
        (1 - (((x ** 2) / r.magnitude) * c(z)))
      end
    end
  end
end
