require "spec_helper"

describe Kepler::Orbit do

  context "when calculating orbital elements from r and v" do
    # values taken from Orbital Mechanics for Engineering Students, Example 4.3
    before(:each) do
      r = Vector[-6045, -3490, 2500]
      v = Vector[-3.457, 6.618, 2.533]
      @orbit = Kepler::Orbit.new(r, v)
    end

    it { expect(@orbit.angular_momentum.round(5)).to eql(Vector[-25385.17, 6669.485, -52070.74]) }
    it { expect(@orbit.radial_velocity.round(5)).to eql(0.55747) }
    it { expect(@orbit.eccentricity.round(5)).to eql(Vector[-0.0916, -0.14221, 0.02644]) }
    it { expect(@orbit.semimajor_axis.round(5)).to eql(8788.09512) }
    it { expect(@orbit.semilatus_rectum.round(5)).to eql(8530.48382) }
    it { expect(@orbit.inclination.round(5)).to eql(153.24923) }
    it { expect(@orbit.node_line.round(5)).to eql(Vector[-6669.485, -25385.17, 0.0]) }
    it { expect(@orbit.right_ascension.round(5)).to eql(255.27929) }
    it { expect(@orbit.argument_of_periapsis.round(5)).to eql(20.06832) }
    it { expect(@orbit.true_anomaly.round(5)).to eql(28.44563) }
    it { expect(@orbit.apoapsis.round(5)).to eql(10292.7255) }
    it { expect(@orbit.periapsis.round(5)).to eql(7283.46473) }
    it { expect(@orbit.period.round(5)).to eql(8198.85762) }
  end

  context "when updating r and v after dt" do
    # values taken from Orbital Mechanics for Engineering Students, Example 3.6
    before(:each) do
      r = Vector[7000, -12124, 0]
      v = Vector[2.6679, 4.6210, 0]
      @orbit = Kepler::Orbit.new(r, v)
    end

    it { expect(@orbit.universal_anomaly(3600).round(5)).to eql(253.53449) }
    it do
      @orbit.update!(3600)
      expect(@orbit.r.round(5)).to eql(Vector[-3297.76863, 7413.39665, 0.0])
      expect(@orbit.v.round(5)).to eql(Vector[-8.2976, -0.96404, -0.0])
    end
  end

  context "set r and v from common params" do
    before(:each) do
      # values taken from Orbital Mechanics for Engineering Students, Example 4.7
      @orbit = Kepler::Orbit.from_params({
        semilatus_rectum: 16056.196688409433,
        eccentricity: 1.4,
        inclination: 30,
        argument_of_periapsis: 60,
        right_ascension: 40,
        true_anomaly: 30
      })
    end

    it do
      expect(@orbit.r.round(5)).to eql(Vector[-4039.89592, 4814.56048, 3628.6247])
      expect(@orbit.v.round(5)).to eql(Vector[-10.38599, -4.77192, 1.74388])
    end
  end

  context "set r and v from apogee, perigee, and inclination" do
    before(:each) do
      # ISS params from http://www.heavens-above.com/orbit.aspx?satid=25544
      @orbit = Kepler::Orbit.from_params({
        apogee: 416,
        perigee: 405,
        inclination: 51.65,
        right_ascension: 304.0847,
        argument_of_periapsis: 117.7713
      })
    end

    it do
      expect(@orbit.r.round(5)).to eql(Vector[1311.56365, 4699.59865, 4701.88142])
      expect(@orbit.v.round(5)).to eql(Vector[-5.64188, 4.37964, -2.80374])
    end
  end

  context "perifocal position from angular momentum, eccentricity, and true_anomaly" do
    # values taken from Orbital Mechanics for Engineering Students, Example 4.7
    it do
      h = 80000
      e = 1.4
      theta = 30
      r_p = Kepler::Orbit.perifocal_position(h, e, theta)

      expect(r_p.round(5)).to eql(Vector[6284.96235, 3628.6247, 0.0])
    end
  end

  context "perifocal velocity from angular momentum, eccentricity, and true_anomaly" do
    # values taken from Orbital Mechanics for Engineering Students, Example 4.7
    it do
      h = 80000
      e = 1.4
      theta = 30
      v_p = Kepler::Orbit.perifocal_velocity(h, e, theta)

      expect(v_p.round(5)).to eql(Vector[-2.49125, 11.29047, 0.0])
    end
  end
end

