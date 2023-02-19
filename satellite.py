from orbital_params import construct_orbital_params
import numpy as np
import datetime
from pymap3d import ecef2eci, geodetic2ecef
import pyquaternion

"""
:param G: float, гравитационная постоянная  Н*м^2*кг^(-2)
:param ME: float, масса Земли, кг
:param MU: float, гравитационный параметр небесного тела
"""
G = 6.6743e-11
ME = 5.9722e24
MU = G * ME

"""
:param JD_START: float, начальный момент времени в юлианских днях
:param UTC_TIME_START, 
:param SECONDS_PER_DAY, 
:param TARGET_POINT,
"""
JD_START = 8084.05644194318
UTC_TIME_START = datetime.datetime(2022, 2, 18, 13, 21, 16, 584)
SECONDS_PER_DAY = 86400
TARGET_POINT = {'lat': 45.920266, 'lon': 63.342286}


class Satellite():
    def __init__(self, start_pos, start_vel):
        self._orb_params = construct_orbital_params(start_pos, start_vel)
        self._coordinates = start_pos
        self._velocity = start_vel
        self._transform_matrix = None
        self._transform_matrix1 = None

    def _count_position(self):
        x, y, z = (self._orb_params._semi_major_axis * (self._transform_matrix * (np.cos(self._orb_params._ecc_anomaly) - self._orb_params._eccentricity) + self._transform_matrix1 *
                                                        np.sqrt(1 - self._orb_params._eccentricity ** 2) * np.sin(self._orb_params._ecc_anomaly)))
        return float(x), float(y), float(z)

    def _count_velocity(self):
        vel_x, vel_y, vel_z = (np.sqrt(MU / self._orb_params._semi_major_axis) / (1 - self._orb_params._eccentricity * np.cos(self._orb_params._ecc_anomaly))) * (- self._transform_matrix) * np.sin(self._orb_params._ecc_anomaly) + self._transform_matrix1 * (
            np.sqrt(1 - self._orb_params._eccentricity ** 2) * np.cos(self._orb_params._ecc_anomaly))

        return float(vel_x), float(vel_y), float(vel_z)

    def count_quanternion(self):
        pass

    def get_params(self, delta_t):
        self._orb_params.update_anomaly(delta_t)

        julian_date = JD_START + delta_t / (SECONDS_PER_DAY)

        x, y, z = self._count_position()
        vel_x, vel_y, vel_z = self._count_velocity()

        return [julian_date, x, y, z, vel_x, vel_y, vel_z]


def construct_satellite(START_POS, START_VEL):
    satellite = Satellite(START_POS, START_VEL)
    ascending_node = satellite._orb_params._ascending_node
    periapsis_arg = satellite._orb_params._periapsis_arg
    inclination = satellite._orb_params._inclination
    satellite._transform_matrix = np.array([[np.cos(ascending_node) * np.cos(periapsis_arg) - np.sin(ascending_node) * np.sin(periapsis_arg) * np.cos(inclination)],
                                            [np.sin(ascending_node) * np.cos(periapsis_arg) + np.cos(ascending_node)
                                             * np.sin(periapsis_arg) * np.cos(inclination)],
                                            [np.sin(periapsis_arg) * np.sin(inclination)]])
    satellite._transform_matrix1 = np.array([[- np.cos(ascending_node) * np.sin(periapsis_arg) - np.sin(ascending_node) * np.cos(periapsis_arg) * np.cos(inclination)],
                                             [-np.sin(ascending_node) * np.sin(periapsis_arg) + np.cos(ascending_node)
                                              * np.cos(periapsis_arg) * np.cos(inclination)],
                                             [np.cos(periapsis_arg) * np.sin(inclination)]])
    return satellite


def convert_geo2eci(target_point, time):
    x, y, z = geodetic2ecef(target_point['lat'], target_point['lon'], 0)
    x, y, z = ecef2eci(x, y, z, time)
    return {'x': x, 'y': y, 'z': z}
