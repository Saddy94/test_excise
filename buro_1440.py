import numpy as np
import kepler
#НУ положения КА и скорости
position  = [4362521.19692133, -2174459.71448059, 4720847.40402189] 
speed = [5356.39915538069, 4741.41348686709, -2761.5472632395]
# G * M
mu_z = 0.000000000066743016 * 5972200000000000000000000  
# Евклидова норма
r = np.linalg.norm(position,ord = 2)
velocity = np.linalg.norm(speed,ord = 2)
#кеплерова энергия
energy = ((velocity**2)/2) - (mu_z / r)
#орбитальный угловой момент
moment = np.linalg.norm(np.cross(position, speed))
#большая полуось
a = - mu_z / (2* energy)
# фокальное расстояние
p = moment**2/(mu_z)
# эксцентриситет
e = np.sqrt((1 + (2 * energy * p)/ (mu_z)))
# период обращения 
T = 2 * np.pi * np.sqrt(a**3/mu_z)
# эксцентрическая аномалия
E = np.arctan2(np.dot(position, speed)/ (np.sqrt(mu_z*a)), 1 - r)
# начальное значение средней аномалии
M = E - np.dot(position, speed)/ (np.sqrt(mu_z*a))
M0 = M - (np.sqrt(mu_z/(a**3)))
# численное решение уравнения Кеплера
eccentric_anomaly, cos_true_anomaly, sin_true_anomaly = kepler.kepler(M0, e)

print(2*np.pi/T, (np.sqrt(mu_z/(a**3))))