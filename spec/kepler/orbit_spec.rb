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
end

