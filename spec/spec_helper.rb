require "rspec"
require "pry-byebug"
require "kepler"

# adding a #round method to Vector so we can set expectations on float components
class Vector
  def round(ndigits = 0)
    map { |e| e.round(ndigits) }
  end
end
