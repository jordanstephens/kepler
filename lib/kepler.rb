require "matrix"

require "kepler/core_ext/numeric"
require "kepler/version"
require "kepler/orbit"

module Kepler
  MU = 398600.0 # km^3s^-2
  EARTH_RADIUS = 6371.0 # km
  I = Vector[1, 0, 0]
  J = Vector[0, 1, 0]
  K = Vector[0, 0, 1]
end
