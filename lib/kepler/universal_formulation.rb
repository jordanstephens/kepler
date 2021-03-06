require "kepler/stumpff"

module Kepler
  # Kepler's equation in terms of the Universal Variable
  module UniversalFormulation
    class << self
      include Stumpff

      def z(x, a)
        (x ** 2) / a
      end

      def f(x, a, r, v, mu, dt)
        z = z(x, a)
        s = s(z)
        c = c(z)

        ((1 - (r.magnitude / a)) * s * (x ** 3)) +
        ((r.inner_product(v) / Math.sqrt(mu)) * c * (x ** 2)) +
        (r.magnitude * x) -
        (Math.sqrt(mu) * dt)
      end

      def dfdt(x, a, r, v, mu)
        z = z(x, a)
        s = s(z)
        c = c(z)

        (c * (x ** 2)) +
        ((r.inner_product(v) / Math.sqrt(mu)) * (1 - (s * z)) * x) +
        (r.magnitude * (1 - (c * z)))
      end

      def d2fdt(x, a, r, v, mu)
        z = z(x, a)
        s = s(z)
        c = c(z)

        ((1 - (r.magnitude / a)) * (1 - (s * z)) * x) +
        ((r.inner_product(v) / Math.sqrt(mu)) * (1 - (c * z)))
      end
    end
  end
end
