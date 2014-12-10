require "kepler/stumpff"

module Kepler
  module UniversalFormulation
    include Stumpff

    def Z(x)
      (x ** 2) / a
    end

    def uf_F(x, dt)
      z = Z(x)
      s = S(z)
      c = C(z)

      ((1 - (r.magnitude / a)) * s * (x ** 3)) +
      ((r.inner_product(v) / Math.sqrt(MU)) * c * (x ** 2)) +
      (r.magnitude * x) -
      (Math.sqrt(MU) * dt)
    end

    def uf_dFdt(x)
      z = Z(x)
      s = S(z)
      c = C(z)

      (c * (x ** 2)) +
      ((r.inner_product(v) / Math.sqrt(MU)) * (1 - (s * z)) * x) +
      (r.magnitude * (1 - (c * z)))
    end

    def uf_d2Fdt(x)
      z = Z(x)
      s = S(z)
      c = C(z)

      ((1 - (r.magnitude / a)) * (1 - (s * z)) * x) +
      ((r.inner_product(v) / Math.sqrt(MU)) * (1 - (c * z)))
    end
  end
end
