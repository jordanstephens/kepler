
module Kepler
  module ParamHelper
    MINIMUM_PARAM_SETS = [
      %i(semimajor_axis eccentricity),
      %i(semilatus_rectum eccentricity),
      %i(apogee perigee)
    ]

    DEFAULT_PARAMS = {
      inclination: 0,
      argument_of_periapsis: 0,
      right_ascension: 0,
      true_anomaly: 0,
      eccentricity: 0
    }

    def expanded_params(params)
      param_set = MINIMUM_PARAM_SETS.find { |rp| rp & params.keys == rp }
      raise ArgumentError, "Invalid parameter set" if param_set.nil?

      params = DEFAULT_PARAMS.merge(params)

      if param_set.include?(:apogee) && param_set.include?(:perigee)
        params[:semimajor_axis] = semimajor_axis_from_apogee_and_perigee(params[:apogee], params[:perigee])
        params[:eccentricity] = eccentricity_from_semimajor_axis_and_perigee(params[:semimajor_axis], params[:perigee])
      elsif param_set.include?(:semilatus_rectum)
        params[:semimajor_axis] = semimajor_axis_from_semilatus_rectum_and_eccentricity(params[:semilatus_rectum], params[:eccentricity])
      end

      unless param_set.include?(:semilatus_rectum)
        params[:semilatus_rectum] = semilatus_rectum_from_semimajor_axis_and_eccentricity(params[:semimajor_axis], params[:eccentricity])
      end

      params[:angular_momentum] = angular_momentum_from_semilatus_rectum(params[:semilatus_rectum])

      params
    end

    def angular_momentum_from_semilatus_rectum(semilatus_rectum)
      Math.sqrt(semilatus_rectum * MU)
    end

    def semilatus_rectum_from_semimajor_axis_and_eccentricity(semimajor_axis, eccentricity)
      semimajor_axis * (1 - (eccentricity ** 2))
    end

    def semimajor_axis_from_semilatus_rectum_and_eccentricity(semilatus_rectum, eccentricity)
      semilatus_rectum / (1 - (eccentricity ** 2))
    end

    def semimajor_axis_from_apogee_and_perigee(apogee, perigee)
      ((EARTH_RADIUS * 2) + apogee + perigee) / 2
    end

    def eccentricity_from_semimajor_axis_and_perigee(semimajor_axis, perigee)
      (semimajor_axis / (EARTH_RADIUS + perigee)) - 1
    end

    # matrix to transform vectors with perifocal basis to vectors with
    # geocentric equatorial basis
    # all args are scalar
    def transform_matrix(argument_of_periapsis, inclination, right_ascension)
      w = argument_of_periapsis.to_rad
      i = inclination.to_rad
      omega = right_ascension.to_rad

      sin_omega = Math.sin(omega)
      cos_omega = Math.cos(omega)
      sin_i = Math.sin(i)
      cos_i = Math.cos(i)
      sin_w = Math.sin(w)
      cos_w = Math.cos(w)

      Matrix[
        [-sin_omega * cos_i * sin_w + (cos_omega * cos_w),
         -sin_omega * cos_i * cos_w - (cos_omega * sin_w),
         sin_omega * sin_i],
        [cos_omega * cos_i * sin_w + (sin_omega * cos_w),
         cos_omega * cos_i * cos_w - (sin_omega * sin_w),
         -cos_omega * sin_i],
        [sin_i * sin_w,
         sin_i * cos_w,
         cos_i]
      ]
    end

    # all args are scalar
    def perifocal_position(angular_momentum, eccentricity, true_anomaly)
      h = angular_momentum
      e = eccentricity
      theta = true_anomaly

      ((h ** 2) / MU) *
      (1 / (1 + (e * Math.cos(theta.to_rad)))) *
      Vector[Math.cos(theta.to_rad), Math.sin(theta.to_rad), 0]
    end

    # all args are scalar
    def perifocal_velocity(angular_momentum, eccentricity, true_anomaly)
      h = angular_momentum
      e = eccentricity
      theta = true_anomaly.to_rad

      (MU / h) *
      Vector[-Math.sin(theta), e + Math.cos(theta), 0]
    end
  end
end
