import numpy as np
import kepler
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pymap3d import ecef2eci
from utilities import llh2xyz

"""
:param G: float, гравитационная постоянная 
:param ME: float, масса Земли
:param MU: float, гравитационный параметр небесного тела
:param TIME: int, интервал времени в секундах, в течение которого необходимо выполнять расчет 
:param position_init: list[float], начальное положение ка 
:param velocity_init: list[float], начальные составляющие вектор скорости ка
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
TIME = 20000

# НУ положения, скорости КА и времени
position_init = [4362521.19692133, -2174459.71448059, 4720847.40402189]
velocity_init = [5356.39915538069, 4741.41348686709, -2761.5472632395]
t0 = 0

# Евклидова норма
r = np.linalg.norm(position_init, ord=2)
velocity_norm = np.linalg.norm(velocity_init, ord=2)
# энергия кеплерова движения
h = ((velocity_norm ** 2) / 2) - (MU / r)
# орбитальный угловой момент
omega = np.cross(position_init, velocity_init)
omega_norm = np.linalg.norm(omega)
# большая полуось
a = - MU / (2 * h)
# фокальное расстояние
p = omega_norm ** 2 / (MU)
# эксцентриситет
e = np.sqrt((1 + (2 * h) * (omega_norm ** 2 / MU ** 2)))
# вычисление угла наклона плоскости орбиты к плоскости экватора
i = np.arccos(omega[2] / omega_norm)
# вычисление долготы восходящего угла 
# w = np.arctan(-(omega[0]/omega[1]))
w = np.arctan2(omega[0], -omega[1])
# Вычисление истинной аномалии
true_anomaly = np.arctan2(np.dot(position_init, velocity_init) / omega_norm,
                          1 - r / p)  # np.arccos(np.dot(ex,position_init) / (r * np.linalg.norm(ex, ord=2)))
# Вычисление аргумента широты
peric = np.arctan2(position_init[2],
                   (position_init[1] * omega[0] - position_init[0] * omega[1]) / omega_norm) - true_anomaly
# Вычисление эксцентрической аномалии
E0 = np.arctan2(np.dot(position_init, velocity_init) / (np.sqrt(MU * a)), 1 - r / a)
# Вычисление средней аномалии на момент времени t0
M0 = E0 - np.dot(position_init, velocity_init) / (np.sqrt(MU * a))
# период обращения 
T = 2 * np.pi * np.sqrt(a ** 3 / MU)

coordinates_X = []
coordinates_Y = []
coordinates_Z = []
velocity_X = []
velocity_Y = []
velocity_Z = []
Ecc = []
# julian_date =
delta_t = 0
# матрицы пересчета координат и скоростей
matrix = np.array([[np.cos(w) * np.cos(peric) - np.sin(w) * np.sin(peric) * np.cos(i)],
                   [np.sin(w) * np.cos(peric) + np.cos(w) * np.sin(peric) * np.cos(i)],
                   [np.sin(peric) * np.sin(i)]])
matrix1 = np.array([[- np.cos(w) * np.sin(peric) - np.sin(w) * np.cos(peric) * np.cos(i)],
                    [-np.sin(w) * np.sin(peric) + np.cos(w) * np.cos(peric) * np.cos(i)],
                    [np.cos(peric) * np.sin(i)]])

Y = 2022
M = 2
D = 18
t = 13 + (21 / 60) + 16.584 / 3600
JD = 367 * Y - int(7 * (Y + int((M + 9) / 12)) / 4) - int(3 * (int((Y + (M - 9) / 7) / 100) + 1) / 4) + int(
    275 * M / 9) + D + 1721028.5 + t / 24
JD2020 = JD - 2451545.0
# расчет параметров траекторного движения на заданном временном инетрвале
with open('result.txt', 'w') as f:
    for i in range(int(TIME)):
        M = M0 + delta_t * np.sqrt(MU / a ** 3)
        E, cos_true_anomaly, sin_true_anomaly = kepler.kepler(M, e)
        delta_t += 1
        X, Y, Z = a * (matrix * (np.cos(E) - e) + matrix1 * (np.sqrt(1 - e ** 2) * np.sin(E)))
        coordinates_X.append(float(X))
        coordinates_Y.append(float(Y))
        coordinates_Z.append(float(Z))
        vel_x, vel_y, vel_z = (np.sqrt(MU / a) / (1 - e * np.cos(E))) * (- matrix) * np.sin(E) + matrix1 * (
                np.sqrt(1 - e ** 2) * np.cos(E))
        velocity_X.append(float(vel_x))
        velocity_Y.append(float(vel_y))
        velocity_Z.append(float(vel_z))
        julian_date = str(JD2020 + delta_t / (3600 * 24))
        f.write(julian_date)
# print(E)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(coordinates_X, coordinates_Y, coordinates_Z)
# ax.plot(vel_x, coordinates_Y, coordinates_Z)
# plt.plot(velocity_Y)
# plt.plot(Ecc)
# plt.show()
