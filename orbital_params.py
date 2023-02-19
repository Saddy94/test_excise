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
    def __init__(self):
        self._semi_major_axis = None
        self._eccentricity = None
        self._inclination = None
        self._ascending_node = None
        self._true_anomaly = None
        self._periapsis_arg = None
        self._mean_anomaly_begin = None
        self._average_distance = None
        self._ecc_anomaly = None
        self._mean_anomaly = None

    def update_anomaly(self, time):
        self._mean_anomaly = self._mean_anomaly_begin + time * self._average_distance
        self._ecc_anomaly, _, _ = kepler.kepler(
            self._mean_anomaly, self._eccentricity)


def construct_orbital_params(start_pos, start_velocity):
    orb_params = OrbitalParams()

    pos_proj = list(start_pos.values())
    vel_proj = list(start_velocity.values())

    pos_norm = np.linalg.norm(pos_proj, ord=2)
    vel_norm = np.linalg.norm(vel_proj, ord=2)

    orb_moment = count_orb_moment(
        pos_proj, vel_proj)
    orb_moment_norm = np.linalg.norm(orb_moment)

    keplerian_energy = count_keplerian_energy(pos_norm, vel_norm)

    focal_distance = count_focal_distance(orb_moment_norm)

    orb_params._semi_major_axis = count_semi_major_axis(keplerian_energy)

    orb_params._eccentricity = count_eccentricity(
        keplerian_energy, orb_moment_norm)

    orb_params._true_anomaly = count_true_anomaly(
        pos_proj, vel_proj, orb_moment_norm, pos_norm, focal_distance)

    orb_params._periapsis_arg = count_periapsis_arg(
        start_pos, orb_moment, orb_moment_norm, orb_params._true_anomaly)

    orb_params._inclination = count_inclination(
        orb_moment, orb_moment_norm)

    orb_params._ascending_node = count_ascending_node(orb_moment)

    ecc_anomaly_begin = count_ecc_anomaly_begin(
        pos_proj, vel_proj, orb_params._semi_major_axis, pos_norm)

    orb_params._mean_anomaly_begin = count_mean_anomaly_begin(
        ecc_anomaly_begin, pos_proj, vel_proj, orb_params._semi_major_axis)

    orb_params._average_distance = count_average_distance(
        orb_params._semi_major_axis)

    orb_params._mean_anomaly = orb_params._mean_anomaly_begin
    orb_params._ecc_anomaly = ecc_anomaly_begin

    return orb_params


def count_keplerian_energy(pos_norm, vel_norm):
    return ((vel_norm ** 2) / 2) - (MU / pos_norm)


def count_orb_moment(pos_proj, vel_proj):
    return np.cross(pos_proj, vel_proj)


def count_eccentricity(keplerian_energy, orb_moment_norm):
    return np.sqrt((1 + (2 * keplerian_energy) * (orb_moment_norm ** 2 / MU ** 2)))


def count_true_anomaly(pos_proj, vel_proj, orb_moment_norm, pos_norm, focal_distance):
    return np.arctan2(np.dot(pos_proj, vel_proj) /
                      orb_moment_norm, 1 - pos_norm / focal_distance)


def count_periapsis_arg(start_pos, orb_moment, orb_moment_norm, true_anomaly):
    return np.arctan2(start_pos['z'], (start_pos['y'] * orb_moment[0] - start_pos['x'] * orb_moment[1]) / orb_moment_norm) - true_anomaly


def count_semi_major_axis(keplerian_energy):
    return - MU / (2 * keplerian_energy)


def count_focal_distance(orb_moment_norm):
    return orb_moment_norm ** 2 / MU


def count_inclination(orb_moment, orb_moment_norm):
    return np.arccos(orb_moment[2] / orb_moment_norm)


def count_ascending_node(orb_moment):
    return np.arctan2(orb_moment[0], - orb_moment[1])


def count_ecc_anomaly_begin(pos_proj, vel_proj, semi_major_axis, pos_norm):
    return np.arctan2(np.dot(pos_proj, vel_proj) / (np.sqrt(MU * semi_major_axis)), 1 - pos_norm / semi_major_axis)


def count_mean_anomaly_begin(ecc_anomaly_begin, pos_proj, vel_proj, semi_major_axis):
    return ecc_anomaly_begin - np.dot(pos_proj, vel_proj) / (np.sqrt(MU * semi_major_axis))


def count_average_distance(semi_major_axis):
    return np.sqrt(MU / semi_major_axis ** 3)
