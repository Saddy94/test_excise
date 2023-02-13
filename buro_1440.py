import numpy as np
import kepler
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#НУ положения, скорости КА и времени
position  = [4362521.19692133, -2174459.71448059, 4720847.40402189] 
speed = [5356.39915538069, 4741.41348686709, -2761.5472632395]
t0 = 0#86400 *  8084.05644194318
# G * M
mu_z = 0.000000000066743016 * 5972200000000000000000000  
# Евклидова норма
r = np.linalg.norm(position,ord = 2)
velocity = np.linalg.norm(speed,ord = 2)
#кеплерова энергия
h = ((velocity**2)/2) - (mu_z / r)
#орбитальный угловой момент
omega = np.cross(position, speed)
omega_norm = np.linalg.norm(omega)
#большая полуось
a = - mu_z / (2* h)
# фокальное расстояние
p = omega_norm**2/(mu_z)
# эксцентриситет
e = np.sqrt((1 + (2 * h ) * (omega_norm**2/mu_z**2)))
# вычисление угла наклона плоскости орбиты к плоскости экватора
i = np.arccos(omega[2]/omega_norm)
# вычисление долготы восходящего угла 
#w = np.arctan(-(omega[0]/omega[1]))
w = np.arctan2(omega[0] , -omega[1])
# вычисление аргумента перицентра
ex = np.cross(speed, omega) / mu_z - position/r
hx = [-omega[1], omega[0], 0]
# peric = np.arccos(np.dot(hx,ex)/(np.linalg.norm(hx, ord=2) * np.linalg.norm(ex, ord=2)))
# Вычисление истинной аномалии
true_anomaly = np.arctan2(np.dot(position, speed)/omega_norm, 1 - r/p)#np.arccos(np.dot(ex,position) / (r * np.linalg.norm(ex, ord=2)))
# Вычисление аргумента широты
peric = np.arctan2(position[2], (position[1] * omega[0] - position[0] * omega[1])/ omega_norm) - true_anomaly
# Вычисление эксцентрической аномалии
E0 = np.arctan2(np.dot(position, speed) / (np.sqrt(mu_z*a)), 1 - r/a)
# Вычисление средней аномалии на момент времени t0
M0 = E0 - np.dot(position, speed)/(np.sqrt(mu_z*a))
# период обращения 
T = 2 * np.pi * np.sqrt(a**3/mu_z)
t = 0*T
delta_t = (t - t0) 
# while delta_t >= T:
#     delta_t -= T
#delta_t = 0  
coordinates_X = []
coordinates_Y = []
coordinates_Z = []
#for i in range((int(T))):

# M = M0 + delta_t * np.sqrt(mu_z/a**3)
# E, cos_true_anomaly, sin_true_anomaly = kepler.kepler(M, e)
#true_anomaly = E + 2 * np.arctan(e * np.sin(E)/ (1 + np.sqrt(1 - e**2) - e * np.cos(E)))


matrix = np.array([[np.cos(w) * np.cos(peric) - np.sin(w) * np.sin(peric) * np.cos(i)],
 				[np.sin(w) * np.cos(peric) + np.cos(w) * np.sin(peric) * np.cos(i)],
 				[np.sin(peric) * np.sin(i)]])
matrix1 = np.array([[- np.cos(w) * np.sin(peric) - np.sin(w) * np.cos(peric) * np.cos(i)],
 				[-np.sin(w) * np.sin(peric) + np.cos(w) * np.cos(peric) * np.cos(i)],
 				[np.cos(peric) * np.sin(i)]])
# delta_t += 1
X,Y,Z = a *( matrix * (np.cos(E0) - e) + matrix1 * (np.sqrt(1 - e**2) * np.sin(E0)))
print(X, Y, Z)



# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# #ax.plot(coordinates_X, coordinates_Y, coordinates_Z)
# plt.plot(coordinates_X, coordinates_Y)
# plt.show()

# Y = 2022
# M = 2
# D = 18
# TIME = (13) + (21/60) + (17.584)/(3600) + 7.7375745507902958333333333333333
# JD = 367* Y - int(7*(Y + int((M + 9)/12))/4) - int(3*(int((Y + (M - 9)/7)/100) + 1)/4) + int(275*M/9) + D + 1721028.5 + TIME/24 
# JD2020 = JD - 2451545.0
# print(JD2020)




