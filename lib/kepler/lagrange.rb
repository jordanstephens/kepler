
module Kepler
  module Lagrange
    include Stumpff

    def f(x, z)
      1 - ((x ** 2) / r.magnitude) * C(z)
    end

    def g(x, z, dt)
      dt - ((1 / Math.sqrt(MU)) * (x ** 3) * S(z))
    end

    def df(x, z, r0)
      (Math.sqrt(MU) / (r.magnitude * r0.magnitude)) * (S(z) * z - 1) * x
    end

    def dg(x, z)
      (1 - (((x ** 2) / r.magnitude) * C(z)))
    end
  end
end
