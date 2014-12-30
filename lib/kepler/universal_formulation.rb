require "kepler/stumpff"

module Kepler
  # Kepler's equation in terms of the Universal Variable
  module UniversalFormulation
    class << self
      include Stumpff

      def z(x, a)
        (x ** 2) / a
      end

      def f(x, a, r, v, dt)
        z = z(x, a)
        s = s(z)
        c = c(z)

        ((1 - (r.magnitude / a)) * s * (x ** 3)) +
        ((r.inner_product(v) / Math.sqrt(MU)) * c * (x ** 2)) +
        (r.magnitude * x) -
        (Math.sqrt(MU) * dt)
      end

      def dfdt(x, a, r, v)
        z = z(x, a)
        s = s(z)
        c = c(z)

        (c * (x ** 2)) +
        ((r.inner_product(v) / Math.sqrt(MU)) * (1 - (s * z)) * x) +
        (r.magnitude * (1 - (c * z)))
      end

      def d2fdt(x, a, r, v)
        z = z(x, a)
        s = s(z)
        c = c(z)

        ((1 - (r.magnitude / a)) * (1 - (s * z)) * x) +
        ((r.inner_product(v) / Math.sqrt(MU)) * (1 - (c * z)))
      end
    end
  end
end
