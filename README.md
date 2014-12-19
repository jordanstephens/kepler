# Kepler

A ruby gem for working with two-body [Keplerian Orbits][0].

## Installation

Release pending, but you can use source it in your Gemfile from Github for now.

**Gemfile**

    gem "kepler", github: "jordanstephens/kepler"

## Usage

### Defining an Orbit

The easiest way to define an orbit is with [common orbital elements][1].

    qzss = Kepler::Orbit.from_params({
      semimajor_axis: 42164, # km
      eccentricity: 0.075,
      inclination: 43, # deg
      right_ascension: 195, # deg
      argument_of_periapsis: 270 # deg
    })

*QZSS orbital parameters from [Wikipedia][2].*

You can also use more colloquial elements like `perigee` and `apogee` instead of `semimajor_axis` and `eccentricity`.

    iss = Kepler::Orbit.from_params({
      apogee: 426.9, # km
      perigee: 416.2, # km
      inclination: 51.65, # deg
      right_ascension: 304.1, # deg
      argument_of_periapsis: 117.8 # deg
    })

*ISS orbital parameters from [Wolfram Alpha][3].*

Or you can use *position* (`r`) and *velocity* (`v`) vectors.

    r = Vector[-6045, -3490, 2500]
    v = Vector[-3.457, 6.618, 2.533]
    orbit = Kepler::Orbit.new(r, v)

### Orbital Attributes

Defining an orbit gives you access to the object's current *position* (`r`) and *velocity* (`v`) vectors along with many other attributes, including:

* `semimajor_axis`
* `semilatus_rectum`
* `apoapsis`
* `periapsis`
* `eccentricity`
* `angular_momentum`
* `radial_velocity`
* `inclination`
* `right_ascension`
* `argument_of_periapsis`
* `true_anomaly`
* `period`

### Updating an Orbit After `dt`

Once you have an initial orbit defined, you can get updated *position* (`r`) and *velocity* (`v`) vectors after a period of time (in seconds) has passed.

    iss.r # => Vector[4427.6294614883145, 662.2433879237291, 5103.355851378229]
    iss.v # => Vector[-3.000020372892516, 6.843005359324254, 1.7091866273043546]

    # Update vectors after 45 minutes have passed
    iss.update!(60 * 45)

	iss.r # => Vector[-4662.9620885090435, -88.56494596532605, -4943.141439770972]
	iss.v # => Vector[2.5073027532818064, -6.876096362825007, -2.248327189819066]

## Caveats

* Only point objects are considered in a two-body environment.
* Perturbations due to atmospheric drag, solar radiation, etc are not considered.
* Nodal precession is not considered.
* Ruby 2.1 or greater is required because STL Vector#cross_product was not defined until Ruby 2.1

## Specs

Specs can be run with

    $ rspec spec

*Most of the specs contain values used in examples found in [Orbital Mechanics for Engineering Students][4].*

[0]: http://en.wikipedia.org/wiki/Kepler_orbit
[1]: http://en.wikipedia.org/wiki/Orbital_elements
[2]: http://en.wikipedia.org/wiki/Quasi-Zenith_Satellite_System
[3]: http://www.wolframalpha.com/input/?i=ISS+orbit
[4]: http://booksite.elsevier.com/9780123747785/
