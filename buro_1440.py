import numpy as np
import kepler
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pymap3d import ecef2eci, geodetic2ecef
from utilities import llh2xyz

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


start_pos = {'x': 4362521.19692133, 'y': -
             2174459.71448059, 'z': 4720847.40402189}
start_velocity = {'vel_x': 5356.39915538069,
                  'vel_y': 4741.41348686709, 'vel_z': -2761.5472632395}
target_point = {'lat': 45.920266, 'lon': 63.342286}
jd_start = 8084.05644194318

# position_norm = np.linalg.norm([*start_pos.values()], ord=2)
# velocity_norm = np.linalg.norm([*start_velocity.values()], ord=2)


class Orbit():
    """
    """
    def __init__(self, start_pos, start_velocity):
        self._start_pos = start_pos
        self._start_velocity = start_velocity
        self._position_norm = np.linalg.norm([*start_pos.values()], ord=2)
        self._velocity_norm = np.linalg.norm([*start_velocity.values()], ord=2)
        self._orbital_momentum = None
        self._keplerian_energy = None
        self._orbital_momentum_norm = None
        self._focal_distance = None
        self._semi_major_axis = None
        self._eccentricity = None
        self._inclination = None
        self._ascending_node = None
        self._true_anomaly = None
        self._periapsis_arg = None
        self._ecc_anomaly = None
        self._mean_anomaly = None
        self._average_distance = None

    def keplerian_energy(self):
        if not self._keplerian_energy:
            self._keplerian_energy = ((self._velocity_norm ** 2) / 2) - \
                (MU / self._position_norm)
        return self._keplerian_energy

    def orbital_momentum(self):
        if not isinstance(self._orbital_momentum, np.ndarray):
            self._orbital_momentum = np.cross(
                list(self._start_pos.values()), list(self._start_velocity.values()))
        return self._orbital_momentum

    def orbital_momentum_norm(self):
        if not self._orbital_momentum_norm:
            self._orbital_momentum_norm = np.linalg.norm(
                self.orbital_momentum())
        return self._orbital_momentum_norm

    def focal_distance(self):
        if not self._focal_distance:
            self._focal_distance = self.orbital_momentum_norm() ** 2 / MU
        return self._focal_distance

    def semi_major_axis(self):
        if not self._semi_major_axis:
            self._semi_major_axis = - MU / (2 * self.keplerian_energy())
        return self._semi_major_axis

    def eccentricity(self):
        if not self._eccentricity:
            self._eccentricity = np.sqrt(
                (1 + (2 * self.keplerian_energy()) * (self.orbital_momentum_norm() ** 2 / MU ** 2)))
        return self._eccentricity

    def inclination(self):
        if not self._inclination:
            self._inclination = np.arccos(self.orbital_momentum()[
                2] / self.orbital_momentum_norm())
        return self._inclination

    def ascending_node(self):
        if not self._ascending_node:
            self._ascending_node = np.arctan2(self.orbital_momentum()[
                0], - self.orbital_momentum()[1])
        return self._ascending_node

    def true_anomaly(self):
        if not self._true_anomaly:
            self._true_anomaly = np.arctan2(np.dot(list(self._start_pos.values()), list(self._start_velocity.values())) / self.orbital_momentum_norm(),
                                            1 - self._position_norm / self.focal_distance())
        return self._true_anomaly

    def periapsis_arg(self):
        if not self._periapsis_arg:
            self._periapsis_arg = np.arctan2(self._start_pos['z'],
                                             (self._start_pos['y'] * self.orbital_momentum()[0] - self._start_pos['x'] * self.orbital_momentum()[1]) / self.orbital_momentum_norm()) - self.true_anomaly()
        return self._periapsis_arg

    def ecc_anomaly(self):
        if not self._ecc_anomaly:
            self._ecc_anomaly = np.arctan2(np.dot(list(self._start_pos.values()), list(self._start_velocity.values())) /
                                           (np.sqrt(MU * self.semi_major_axis())), 1 - self._position_norm / self.semi_major_axis())
        return self._ecc_anomaly

    def mean_anomaly(self):
        if not self._mean_anomaly:
            self._mean_anomaly = self.ecc_anomaly() - np.dot(list(self._start_pos.values()),
                                                             list(self._start_velocity.values()) / (np.sqrt(MU * self.semi_major_axis())))
        return self._mean_anomaly

    def average_distance(self):
        if not self._average_distance:
            self._average_distance = np.sqrt(MU / self.semi_major_axis() ** 3)
        return self._average_distance


coordinates_X = []
coordinates_Y = []
coordinates_Z = []
velocity_X = []
velocity_Y = []
velocity_Z = []
Ecc = []
orbit = Orbit(start_pos, start_velocity)
# матрицы пересчета координат и скоростей
periapsis_arg = orbit.periapsis_arg()
inclination = orbit.inclination()
ascending_node = orbit.ascending_node()
eccentricity = orbit.eccentricity()
semi_major_axis = orbit.semi_major_axis()
true_anomaly = orbit.true_anomaly()
focal_distance = orbit.focal_distance()
# print(periapsis_arg)
print(inclination)
print(ascending_node)
print(true_anomaly)
print(focal_distance)
matrix = np.array([[np.cos(ascending_node) * np.cos(periapsis_arg) - np.sin(ascending_node) * np.sin(periapsis_arg) * np.cos(inclination)],
                   [np.sin(ascending_node) * np.cos(periapsis_arg) + np.cos(ascending_node)
                    * np.sin(periapsis_arg) * np.cos(inclination)],
                   [np.sin(periapsis_arg) * np.sin(inclination)]])
matrix1 = np.array([[- np.cos(ascending_node) * np.sin(periapsis_arg) - np.sin(ascending_node) * np.cos(periapsis_arg) * np.cos(inclination)],
                    [-np.sin(ascending_node) * np.sin(periapsis_arg) + np.cos(ascending_node)
                     * np.cos(periapsis_arg) * np.cos(inclination)],
                    [np.cos(periapsis_arg) * np.sin(inclination)]])

Y = 2022
M = 2
D = 18
t = 13 + (21 / 60) + 16.584 / 3600
JD = 367 * Y - int(7 * (Y + int((M + 9) / 12)) / 4) - int(3 * (int((Y + (M - 9) / 7) / 100) + 1) / 4) + int(
    275 * M / 9) + D + 1721028.5 + t / 24
JD2020 = JD - 2451545.0
# расчет параметров траекторного движения на заданном временном инетрвале


#  TODO: разобраться с шагом, чтобы можно было разные шаги делать, а не только delta_t=1
def process(time_begin, time_end, delta_t=1):
    with open('result.txt', 'w') as f:
        for i in range(int(time_end)):
            M = orbit.mean_anomaly() + 0 * orbit.average_distance()
            E, cos_true_anomaly, sin_true_anomaly = kepler.kepler(
                M, orbit.eccentricity())
            X, Y, Z = semi_major_axis * (matrix * (np.cos(E) - eccentricity) + matrix1 *
                                         (np.sqrt(1 - eccentricity ** 2) * np.sin(E)))

            # vel_x, vel_y, vel_z = (np.sqrt(MU / a) / (1 - e * np.cos(E))) * (- matrix) * np.sin(E) + matrix1 * (
           #    np.sqrt(1 - e ** 2) * np.cos(E))
            # velocity_X.append(float(vel_x))
            # velocity_Y.append(float(vel_y))
            # velocity_Z.append(float(vel_z))
            coordinates_X.append(float(X))
            coordinates_Y.append(float(Y))
            coordinates_Z.append(float(Z))

            julian_date = str(JD2020 + delta_t / (3600 * 24))
            f.write(julian_date)
        print(X, Y, Z)


if __name__ == "__main__":
    process(TIME_BEGIN, TIME_END)

# # print(E)
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.plot(coordinates_X, coordinates_Y, coordinates_Z)
# ax.plot(vel_x, coordinates_Y, coordinates_Z)
# plt.plot(velocity_Y)
# plt.plot(Ecc)
# plt.show()
