require "spec_helper"

describe Kepler::Orbit do

  context "when calculating orbital elements from r and v" do
    # values taken from Orbital Mechanics for Engineering Students, Example 4.3
    before(:each) do
      r = Vector[-6045, -3490, 2500]
      v = Vector[-3.457, 6.618, 2.533]
      @orbit = Kepler::Orbit.new(r, v)
    end

    it { expect(@orbit.angular_momentum).to eql(Vector[-25385.17, 6669.484999999999, -52070.740000000005]) }
    it { expect(@orbit.radial_velocity).to eql(0.5574679274498466) }
    it { expect(@orbit.eccentricity).to eql(Vector[-0.09160485604616679, -0.14220737156769422, 0.02644392824064549]) }
    it { expect(@orbit.semimajor_axis).to eql(8788.095117377656) }
    it { expect(@orbit.semilatus_rectum).to eql(8530.483818970712) }
    it { expect(@orbit.inclination).to eql(153.2492285182475) }
    it { expect(@orbit.node_line).to eql(Vector[-6669.484999999999, -25385.17, 0.0]) }
    it { expect(@orbit.right_ascension).to eql(255.27928533439615) }
    it { expect(@orbit.argument_of_periapsis).to eql(20.068316650582467) }
    it { expect(@orbit.true_anomaly).to eql(28.445628306614967) }
    it { expect(@orbit.apoapsis).to eql(10292.725501794834) }
    it { expect(@orbit.periapsis).to eql(7283.464732960478) }
    it { expect(@orbit.period).to eql(8198.857616829207) }
  end

  context "when updating r and v after dt" do
    # values taken from Orbital Mechanics for Engineering Students, Example 3.6
    before(:each) do
      r = Vector[7000, -12124, 0]
      v = Vector[2.6679, 4.6210, 0]
      @orbit = Kepler::Orbit.new(r, v)
    end

    it { expect(@orbit.universal_anomaly(3600)).to eql(253.53449076412872) }
    it do
      @orbit.update!(3600)
      expect(@orbit.r).to eql(Vector[-3297.7686251992823, 7413.396645787405, 0.0])
      expect(@orbit.v).to eql(Vector[-8.297603024266525, -0.9640449446737689, -0.0])
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
      expect(@orbit.r).to eql(Vector[-4039.895923201738, 4814.560480182377, 3628.6247021718837])
      expect(@orbit.v).to eql(Vector[-10.385987618194685, -4.771921637340853, 1.7438750000000005])
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
      expect(@orbit.r).to eql(Vector[1311.5636463724593, 4699.598648496937, 4701.881415410206])
      expect(@orbit.v).to eql(Vector[-5.641883624845517, 4.379639136504691, -2.803740788571159])
    end
  end

  context "perifocal position from angular momentum, eccentricity, and true_anomaly" do
    # values taken from Orbital Mechanics for Engineering Students, Example 4.7
    it do
      h = 80000
      e = 1.4
      theta = 30
      r_p = Kepler::Orbit.perifocal_position(h, e, theta)

      expect(r_p).to eql(Vector[6284.962345761189, 3628.6247021718837, 0.0])
    end
  end

  context "perifocal velocity from angular momentum, eccentricity, and true_anomaly" do
    # values taken from Orbital Mechanics for Engineering Students, Example 4.7
    it do
      h = 80000
      e = 1.4
      theta = 30
      v_p = Kepler::Orbit.perifocal_velocity(h, e, theta)

      expect(v_p).to eql(Vector[-2.4912499999999995, 11.290471574355966, 0.0])
    end
  end
end

