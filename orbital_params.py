import numpy as np
import kepler

"""
:param G: float, гравитационная постоянная  Н*м^2*кг^(-2)
:param ME: float, масса Земли, кг
:param MU: float, гравитационный параметр небесного тела
"""
G = 6.6743e-11
ME = 5.9722e24
MU = G * ME

"""
:param e: float, эксцентриситет орбиты
:param h: float, кеплерова энергия движения
:param a: float, большая полуось орбиты
:param p: float, фокальное расстояние
:param i: float, угол наклонения орбиты к плоскости экватора
:param true_anomaly: float, истинная аномалия
:param E0: float, эксцентрическая аномалия на эпоху t0
:param M0: float, средняя аномалия на эпоху t0
:param T: float, период обращения спутника по орбите
:param w
"""
class OrbitalParams():
    def __init__(self, start_pos, start_velocity):
        self._start_pos = start_pos
        self._start_velocity = start_velocity
        self._position_norm = None
        self._velocity_norm = None
        self._orbital_momentum = None
        self._orbital_momentum_norm = None
        self._keplerian_energy = None
        self._focal_distance = None
        self._semi_major_axis = None
        self._eccentricity = None
        self._inclination = None
        self._ascending_node = None
        self._true_anomaly = None
        self._periapsis_arg = None
        self._ecc_anomaly_begin = None
        self._mean_anomaly_begin = None
        self._average_distance = None
        self._ecc_anomaly = None
        self._mean_anomaly = None

    def get_orbital_momentum(self):
        return self._orbital_momentum

    def get_orbital_momentum_norm(self):
        return self._orbital_momentum_norm

    def get_keplerian_energy(self):
        return self._keplerian_energy

    def get_focal_distance(self):
        return self._focal_distance

    def get_semi_major_axis(self):
        return self._semi_major_axis

    def get_eccentricity(self):
        return self._eccentricity

    def get_inclination(self):
        return self._inclination

    def get_ascending_node(self):
        return self._ascending_node

    def get_true_anomaly(self):
        return self._true_anomaly

    def get_periapsis_arg(self):
        return self._periapsis_arg

    def get_average_distance(self):
        return self._average_distance

    def update_anomaly(self, time):
        self._mean_anomaly = self._mean_anomaly_begin + time * self._average_distance
        self._ecc_anomaly, _, _ = kepler.kepler(
            self._mean_anomaly, self._eccentricity)

    def get_mean_anomaly(self):
        return self._mean_anomaly

    def get_ecc_anomaly(self):
        return self._ecc_anomaly


def construct_orbital_params(start_pos, start_velocity):
    orbital_params = OrbitalParams(start_pos, start_velocity)

    position_projections = list(orbital_params._start_pos.values())
    velocity_projections = list(orbital_params._start_velocity.values())

    orbital_params._position_norm = np.linalg.norm(position_projections, ord=2)

    orbital_params._velocity_norm = np.linalg.norm(velocity_projections, ord=2)

    orbital_params._orbital_momentum = np.cross(
        position_projections, velocity_projections)

    orbital_params._keplerian_energy = (
        (orbital_params._velocity_norm ** 2) / 2) - (MU / orbital_params._position_norm)

    orbital_params._orbital_momentum_norm = np.linalg.norm(
        orbital_params._orbital_momentum)

    orbital_params._focal_distance = orbital_params._orbital_momentum_norm ** 2 / MU

    orbital_params._semi_major_axis = - MU / \
        (2 * orbital_params._keplerian_energy)

    orbital_params._eccentricity = np.sqrt(
        (1 + (2 * orbital_params._keplerian_energy) * (orbital_params._orbital_momentum_norm ** 2 / MU ** 2)))

    orbital_params._inclination = np.arccos(
        orbital_params._orbital_momentum[2] / orbital_params._orbital_momentum_norm)

    orbital_params._ascending_node = np.arctan2(
        orbital_params._orbital_momentum[0], - orbital_params._orbital_momentum[1])

    orbital_params._true_anomaly = np.arctan2(np.dot(position_projections, velocity_projections) /
                                              orbital_params._orbital_momentum_norm, 1 - orbital_params._position_norm / orbital_params._focal_distance)

    orbital_params._periapsis_arg = np.arctan2(orbital_params._start_pos['z'], (orbital_params._start_pos['y'] * orbital_params._orbital_momentum[0] -
                                                                                orbital_params._start_pos['x'] * orbital_params._orbital_momentum[1]) / orbital_params._orbital_momentum_norm) - orbital_params._true_anomaly

    orbital_params._ecc_anomaly_begin = np.arctan2(np.dot(position_projections, velocity_projections) / (np.sqrt(
        MU * orbital_params._semi_major_axis)), 1 - orbital_params._position_norm / orbital_params._semi_major_axis)

    orbital_params._mean_anomaly_begin = orbital_params._ecc_anomaly_begin - \
        np.dot(position_projections, velocity_projections) / \
        (np.sqrt(MU * orbital_params._semi_major_axis))

    orbital_params._average_distance = np.sqrt(
        MU / orbital_params._semi_major_axis ** 3)

    orbital_params._mean_anomaly = orbital_params._mean_anomaly_begin

    orbital_params._ecc_anomaly = orbital_params._ecc_anomaly_begin

    return orbital_params