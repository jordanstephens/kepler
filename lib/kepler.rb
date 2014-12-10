require "matrix"

require "kepler/version"
require "kepler/orbit"

module Kepler
  MU = 398600.0 # km^3s^-2
  I = Vector[1, 0, 0]
  J = Vector[0, 1, 0]
  K = Vector[0, 0, 1]
end
