module Kepler
  # Laguerre's method for finding roots of a polynomial using its
  # first and second derivatives and an initial guess
  module Laguerre
    def self.solve(guess, pF, pdF, pd2F)

      x = guess
      n = 0

      begin
        f = pF.call(x)
        f_prime = pdF.call(x)
        f_prime_prime = pd2F.call(x)

        delta = 2 * Math.sqrt(
          (4 * (f_prime ** 2)) -
          (5 * f * f_prime_prime)
        )

        dx = (5 * f) / (
          f_prime + ((f_prime.abs / f_prime) * delta)
        )

        x = x - dx
        n += 1
      end until dx == 0 || n > 10

      x
    end
  end
end
