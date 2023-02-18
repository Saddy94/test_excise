import numpy as np
import kepler
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pymap3d import ecef2eci, geodetic2ecef
from utilities import llh2xyz
import datetime
"""
:param G: float, гравитационная постоянная  Н*м^2*кг^(-2)
:param ME: float, масса Земли, кг
:param MU: float, гравитационный параметр небесного тела

:param TIME_BEGIN: int, начальный момент времени в секундах
:param TIME_END: int, конечный момент времени в секундах
:param start_pos: list[float], начальное положение ка 
:param start_velocity: list[float], начальные составляющие вектор скорости ка
:param jd_start: float, начальный момент времени в юлианских днях

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
G = 6.6743e-11
ME = 5.9722e24
MU = G * ME

TIME_BEGIN = 0
TIME_END = 20000

START_POS = {'x': 4362521.19692133, 'y': -
             2174459.71448059, 'z': 4720847.40402189}
START_VEL = {'vel_x': 5356.39915538069,
             'vel_y': 4741.41348686709, 'vel_z': -2761.5472632395}
#time = datetime.datetime(2022, 2, 18, 13, 21, 16, 584)
TARGET_POINT = {'lat': 45.920266, 'lon': 63.342286}
JD_START = 8084.05644194318


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

    orbital_params._keplerian_energy = ((orbital_params._velocity_norm ** 2) / 2) - \
        (MU / orbital_params._position_norm)

    orbital_params._orbital_momentum_norm = np.linalg.norm(
        orbital_params._orbital_momentum)

    orbital_params._focal_distance = orbital_params._orbital_momentum_norm ** 2 / MU

    orbital_params._semi_major_axis = - MU / \
        (2 * orbital_params._keplerian_energy)

    orbital_params._eccentricity = np.sqrt(
        (1 + (2 * orbital_params._keplerian_energy) * (orbital_params._orbital_momentum_norm ** 2 / MU ** 2)))

    orbital_params._inclination = np.arccos(orbital_params._orbital_momentum[
        2] / orbital_params._orbital_momentum_norm)

    orbital_params._ascending_node = np.arctan2(orbital_params._orbital_momentum[
        0], - orbital_params._orbital_momentum[1])

    orbital_params._true_anomaly = np.arctan2(np.dot(position_projections, velocity_projections) / orbital_params._orbital_momentum_norm,
                                              1 - orbital_params._position_norm / orbital_params._focal_distance)

    orbital_params._periapsis_arg = np.arctan2(orbital_params._start_pos['z'],
                                               (orbital_params._start_pos['y'] * orbital_params._orbital_momentum[0] - orbital_params._start_pos['x'] * orbital_params._orbital_momentum[1]) / orbital_params._orbital_momentum_norm) - orbital_params._true_anomaly

    orbital_params._ecc_anomaly_begin = np.arctan2(np.dot(position_projections, velocity_projections) /
                                                   (np.sqrt(MU * orbital_params._semi_major_axis)), 1 - orbital_params._position_norm / orbital_params._semi_major_axis)

    orbital_params._mean_anomaly_begin = orbital_params._ecc_anomaly_begin - np.dot(position_projections,
                                                                                    velocity_projections) / (np.sqrt(MU * orbital_params._semi_major_axis))

    orbital_params._average_distance = np.sqrt(
        MU / orbital_params._semi_major_axis ** 3)

    orbital_params._mean_anomaly = orbital_params._mean_anomaly_begin
    orbital_params._ecc_anomaly = orbital_params._ecc_anomaly_begin

    return orbital_params


class Satellite():
    def __init__(self, start_pos, start_vel):
        self._orbital_params = construct_orbital_params(start_pos, start_vel)
        self._coordinates = start_pos
        self._velocity = start_vel
        self._transform_matrix = None
        self._transform_matrix1 = None

    def get_params(self, delta_t):
        self._orbital_params.update_anomaly(delta_t)

        semi_major_axis = self._orbital_params.get_semi_major_axis()
        ecc_anomaly = self._orbital_params.get_ecc_anomaly()
        eccentricity = self._orbital_params.get_eccentricity()

        x_eci, y_eci, z_eci = semi_major_axis * (self._transform_matrix * (np.cos(ecc_anomaly) - eccentricity) + self._transform_matrix1 *
                                                 np.sqrt(1 - eccentricity ** 2) * np.sin(ecc_anomaly))

        vel_x, vel_y, vel_z = (np.sqrt(MU / semi_major_axis) / (1 - eccentricity * np.cos(ecc_anomaly))) * (- self._transform_matrix) * np.sin(ecc_anomaly) + self._transform_matrix1 * (
            np.sqrt(1 - eccentricity ** 2) * np.cos(ecc_anomaly))
            
        julian_date = JD_START + delta_t / (3600 * 24)

        return [julian_date, float(x_eci), float(y_eci), float(z_eci), float(vel_x), float(vel_y), float(vel_z)]


def construct_satellite(START_POS, START_VEL):
    satellite = Satellite(START_POS, START_VEL)
    ascending_node = satellite._orbital_params.get_ascending_node()
    periapsis_arg = satellite._orbital_params.get_periapsis_arg()
    inclination = satellite._orbital_params.get_inclination()
    satellite._transform_matrix = np.array([[np.cos(ascending_node) * np.cos(periapsis_arg) - np.sin(ascending_node) * np.sin(periapsis_arg) * np.cos(inclination)],
                                            [np.sin(ascending_node) * np.cos(periapsis_arg) + np.cos(ascending_node)
                                             * np.sin(periapsis_arg) * np.cos(inclination)],
                                            [np.sin(periapsis_arg) * np.sin(inclination)]])
    satellite._transform_matrix1 = np.array([[- np.cos(ascending_node) * np.sin(periapsis_arg) - np.sin(ascending_node) * np.cos(periapsis_arg) * np.cos(inclination)],
                                             [-np.sin(ascending_node) * np.sin(periapsis_arg) + np.cos(ascending_node)
                                            * np.cos(periapsis_arg) * np.cos(inclination)],
                                             [np.cos(periapsis_arg) * np.sin(inclination)]])
    return satellite


def convert_data_to_str(data):
    result = ' '.join(str(i) for i in data) + '\n'
    return result

# coordinates_X = []
# coordinates_Y = []
# coordinates_Z = []
# velocity_X = []
# velocity_Y = []
# velocity_Z = []
# Ecc = []
# TODO: разобраться с шагом, чтобы можно было разные шаги делать, а не только delta_t=1
def process(time_begin, time_end, delta_t=1):
    satellite = construct_satellite(START_POS, START_VEL)
    with open('result.txt', 'w') as f:
        for i in range(int(time_end)):
            out_data = satellite.get_params(i)
            out_data = convert_data_to_str(out_data)
            f.writelines(out_data)
            # velocity_X.append(float(vel_x))
            # velocity_Y.append(float(vel_y))
            # velocity_Z.append(float(vel_z))
            # coordinates_X.append(float(x_eci))
            # coordinates_Y.append(float(y_eci))
            # coordinates_Z.append(float(z_eci))

if __name__ == "__main__":
    process(TIME_BEGIN, TIME_END)

    # print(E)
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
#     ax.plot(coordinates_X, coordinates_Y, coordinates_Z)
#     #ax.plot(vel_x, coordinates_Y, coordinates_Z)
#    # plt.plot((coordinates_X, coordinates_Y)
#    # plt.plot(Ecc)
#     plt.show()
