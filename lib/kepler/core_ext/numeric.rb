class Numeric
  def to_deg
    self * 180 / Math::PI
  end

  def to_rad
    self * Math::PI / 180
  end
end
