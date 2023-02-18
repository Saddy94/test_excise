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
TIME_END = 3


START_POS = {'x': 4362521.19692133, 'y': -
             2174459.71448059, 'z': 4720847.40402189}
START_VEL = {'vel_x': 5356.39915538069,
                  'vel_y': 4741.41348686709, 'vel_z': -2761.5472632395}
TARGET_POINT = {'lat': 45.920266, 'lon': 63.342286}
JD_START = 8084.05644194318


class Orbit():
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
        self._ecc_anomaly = None
        self._mean_anomaly = None
        self._average_distance = None

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

    def get_ecc_anomaly(self):
        return self._ecc_anomaly

    def get_mean_anomaly(self):
        return self._mean_anomaly

    def get_average_distance(self):
        return self._average_distance


def construct_orbit(start_pos, start_velocity):
    orbit = Orbit(start_pos, start_velocity)

    position_projections = list(orbit._start_pos.values())
    velocity_projections = list(orbit._start_velocity.values())

    orbit._position_norm = np.linalg.norm(position_projections, ord=2)

    orbit._velocity_norm = np.linalg.norm(velocity_projections, ord=2)

    orbit._orbital_momentum = np.cross(
        position_projections, velocity_projections)

    orbit._keplerian_energy = ((orbit._velocity_norm ** 2) / 2) - \
        (MU / orbit._position_norm)

    orbit._orbital_momentum_norm = np.linalg.norm(
        orbit._orbital_momentum)

    orbit._focal_distance = orbit._orbital_momentum_norm ** 2 / MU

    orbit._semi_major_axis = - MU / (2 * orbit._keplerian_energy)

    orbit._eccentricity = np.sqrt(
        (1 + (2 * orbit._keplerian_energy) * (orbit._orbital_momentum_norm ** 2 / MU ** 2)))

    orbit._inclination = np.arccos(orbit._orbital_momentum[
        2] / orbit._orbital_momentum_norm)

    orbit._ascending_node = np.arctan2(orbit._orbital_momentum[
        0], - orbit._orbital_momentum[1])

    orbit._true_anomaly = np.arctan2(np.dot(position_projections, velocity_projections) / orbit._orbital_momentum_norm,
                                     1 - orbit._position_norm / orbit._focal_distance)

    orbit._periapsis_arg = np.arctan2(orbit._start_pos['z'],
                                      (orbit._start_pos['y'] * orbit._orbital_momentum[0] - orbit._start_pos['x'] * orbit._orbital_momentum[1]) / orbit._orbital_momentum_norm) - orbit._true_anomaly

    orbit._ecc_anomaly = np.arctan2(np.dot(position_projections, velocity_projections) /
                                    (np.sqrt(MU * orbit._semi_major_axis)), 1 - orbit._position_norm / orbit._semi_major_axis)

    orbit._mean_anomaly = orbit._ecc_anomaly - np.dot(position_projections,
                                                      velocity_projections) / (np.sqrt(MU * orbit._semi_major_axis))

    orbit._average_distance = np.sqrt(MU / orbit._semi_major_axis ** 3)

    return orbit


coordinates_X = []
coordinates_Y = []
coordinates_Z = []
velocity_X = []
velocity_Y = []
velocity_Z = []
Ecc = []
orbit = construct_orbit(START_POS, START_VEL)
periapsis_arg = orbit.get_periapsis_arg()
inclination = orbit.get_inclination()
ascending_node = orbit.get_ascending_node()
eccentricity = orbit.get_eccentricity()
semi_major_axis = orbit.get_semi_major_axis()
true_anomaly = orbit.get_true_anomaly()
focal_distance = orbit.get_focal_distance()
keplerian_energy = orbit.get_keplerian_energy()
print(inclination)
print(ascending_node)
print(true_anomaly)
print(semi_major_axis)
print(keplerian_energy)
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
			M = orbit.get_mean_anomaly() + i * orbit.get_average_distance()
			print(M)
			E, cos_true_anomaly, sin_true_anomaly = kepler.kepler(
				M, orbit.get_eccentricity())
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
		#print(len(coordinates_X))#, Y, Z)
		print(coordinates_X[2], coordinates_Y[2], coordinates_Z[2])#, Y, Z)

if __name__ == "__main__":
    process(TIME_BEGIN, TIME_END)

# # print(E)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(coordinates_X, coordinates_Y, coordinates_Z)
#ax.plot(vel_x, coordinates_Y, coordinates_Z)
#plt.plot(velocity_Y)
#plt.plot(Ecc)
plt.show()