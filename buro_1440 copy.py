import numpy as np
import kepler
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#НУ положения, скорости КА и времени
position  = [4362521.19692133, 2174459.71448059, -4720847.40402189] 
speed = [5356.39915538069, 4741.41348686709, -2761.5472632395]
#position  = [2937656.611, 14432705.729, -20836304.223] 
#speed = [-2408.799, 2723.781, 1545.981]
t0 = 0#86400 *  8084.05644194318
# G * M
mu_z = 0.000000000066743016 * 5972200000000000000000000  
# Евклидова норма
r = np.linalg.norm(position,ord = 2)
velocity = np.linalg.norm(speed,ord = 2)
c1 = position[1] * speed[2] - speed[1] * position[2]
c2 = speed[0] * position[2] - position[0] * speed[2] 
c3 = position[0] * speed[1] - speed[0] * position[1]
c = [c1, c2, c3]
# Евклидова норма
c_norm = np.linalg.norm(c, ord = 2)
h = velocity**2 - 2 * mu_z/r
p = c_norm**2 / mu_z
a = - mu_z / h
i = np.arccos(c3/c_norm)
omega = np.arccos(- c2 / (c_norm * np.sin(i)))

D = position[0] * speed[0] + position[1] * speed[1] + position[2] * speed[2]
D_der = mu_z / r + h

f1 = D_der * position[0] - D * speed[0]
f2 = D_der * position[1] - D * speed[1]
f3 = D_der * position[2] - D * speed[2]

f = [f1, f2, f3]
f_norm = np.linalg.norm(f, ord = 2)

e = f_norm / mu_z
per = np.arccos((f1 * np.cos(omega) + f2 * np.sin(omega)) / f_norm)

true_anomaly = np.arccos((p - r)/(e*r))

E = np.arccos((e + np.cos(true_anomaly))/ (1 + e * np.cos(true_anomaly)))

M = E - e * np.sin(E)

u = per + true_anomaly

r_ = a * ((1 - e**2)/ (1 + e * np.cos(true_anomaly))) 

matrix = np.array([[np.cos(omega) * np.cos(u) - np.sin(omega) * np.sin(u) * np.cos(i)],
				[np.sin(omega) * np.cos(u) + np.cos(omega) * np.sin(u) * np.cos(i)],
				[np.sin(u) * np.sin(i)]])

X,Y,Z = r_ * matrix 
# coordinates_X.append(float(X))
# coordinates_Y.append(float(Y))
# coordinates_Z.append(float(Z))
#print(X, Y, Z)
print(i, omega, a)


