import numpy as np
import kepler

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
e = np.sqrt((1 + (2 * h * p)/ (mu_z)))
# вычисление угла наклона плоскости орбиты к плоскости экватора
i = np.arccos(omega[2]/omega_norm)
# вычисление долготы восходящего угла 
w = np.arccos(- omega[1]/ (omega_norm))
# вычисление аргумента перицентра
ex = np.cross(speed, omega) / mu_z - position/r
hx = [-omega[1], omega[0], 0]
peric = np.arccos(np.dot(hx,ex)/(np.linalg.norm(hx, ord=2) * np.linalg.norm(ex, ord=2)))
# Вычисление истинной аномалии
true_anomaly = np.arccos(np.dot(ex,position) / (r * np.linalg.norm(ex, ord=2)))
# Вычисление эксцентрической аномалии
E = 2 * np.arctan(np.tan(true_anomaly/2)/(np.sqrt((1+e) / (1-e))))
# Вычисление средней аномалии на момент времени t0
M0 = E - e * np.sin(E)
# период обращения 
T = 2 * np.pi * np.sqrt(a**3/mu_z)
t = 5*T
delta_t = (t - t0) 
if delta_t >= T:
    delta_t -= int(delta_t/T)
M = M0 + T * np.sqrt(mu_z/a**3)

#time = 8084.05644194318 + int(t/3600) + int()
print(a,e,i,w,peric, true_anomaly, E, M0, M, T)

Y = 2022
M = 2
D = 18
TIME = (13) + (21/60) + (17.584)/(3600) + 7.7375745507902958333333333333333
JD = 367* Y - int(7*(Y + int((M + 9)/12))/4) - int(3*(int((Y + (M - 9)/7)/100) + 1)/4) + int(275*M/9) + D + 1721028.5 + TIME/24 
JD2020 = JD - 2451545.0
print(JD2020)




