module Kepler
  module Stumpff
    def C(z)
      if z > 0
        (1 - Math.cos(Math.sqrt(z))) / z
      elsif z < 0
        (Math.cosh(Math.sqrt(-z)) - 1) / (-z)
      else
        0.5
      end
    end

    def S(z)
      if z > 0
        (Math.sqrt(z) - Math.sin(Math.sqrt(z))) / (Math.sqrt(z) ** 3)
      elsif z < 0
        (Math.sinh(Math.sqrt(-z)) - Math.sqrt(-z))/(Math.sqrt(-z) ** 3)
      else
        1 / 6
      end
    end
  end
end
