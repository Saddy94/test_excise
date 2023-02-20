from orbital_params import construct_orbital_params
import numpy as np
import pyquaternion
from physical_const import *


# начальный момент времени в юлианских днях
JD_START = 8084.05644194318
# количество секунд в сутках
SECONDS_PER_DAY = 86400


def count_position(semi_major_axis, ecc_anomaly, eccentricity, transf_mat):
    """
    расчет положения спутника
    Returns:
        pos (dict): 
            x (float): x - компонента радиус-вектора спутника в осях ИСК
            y (float): y - компонента радиус-вектора спутника в осях ИСК
            z (float): z - компонента радиус-вектора спутника в осях ИСК
    """
    x, y, z = (semi_major_axis * (transf_mat[0]
        * (np.cos(ecc_anomaly) - eccentricity) + transf_mat[1] *
        np.sqrt(1 - eccentricity ** 2) * np.sin(ecc_anomaly)))
    pos = {'x': float(x), 'y': float(y), 'z': float(z)}

    return pos


def count_velocity(semi_major_axis, ecc_anomaly, eccentricity, transf_mat):
    """
    расчет вектора скорости спутника 
    Returns:
        vel (dict): 
            vel_x (float): x - компонента скорости спутника в осях ИСК
            vel_y (float): y - компонента скорости спутника в осях ИСК
            vel_z (float): z - компонента скорости спутника в осях ИСК
    """
    vel_x, vel_y, vel_z = (np.sqrt(MU / semi_major_axis)
        / (1 - eccentricity * np.cos(ecc_anomaly))) * ((- transf_mat[0]) *
        np.sin(ecc_anomaly) + transf_mat[1] * (np.sqrt(1 - eccentricity ** 2)
        * np.cos(ecc_anomaly)))
    vel = {'vel_x': float(vel_x), 'vel_y': float(vel_y), 'vel_z': float(vel_z)}

    return vel


def count_quaternion(true_anomaly, periapsis_arg, inclination, ascending_node,
    pos, target_point):
    """
    расчет квантерниона целевой ориентации спутника в связанной СК
    
    алгоритм выполняется в два этапа: сначала вычисляется 
    матрица, определяющая связь ИСК и орбитальной СК, после чего 
    определяется квантернион наведения оси Y ССК на целевую точку.
    Квантернион целевой ориентации определяется произведением 2 
    вышеназванных квантернионовю
    Returns:
        target_q (pyquaternion.Quaternion): квантернион целевой ориентации КА
    """
    lat_arg = true_anomaly + periapsis_arg

    matrix11 = np.cos(ascending_node) * np.cos(lat_arg) - \
        np.sin(ascending_node) * np.sin(lat_arg) * np.cos(inclination)
    matrix21 = np.sin(ascending_node) * np.cos(lat_arg) + \
        np.cos(ascending_node) * np.sin(lat_arg) * np.cos(inclination)
    matrix31 = np.sin(lat_arg) * np.sin(inclination)

    matrix12 = - np.cos(ascending_node) * np.sin(lat_arg) - \
        np.sin(ascending_node) * np.cos(lat_arg) * np.cos(inclination)
    matrix22 = - np.sin(ascending_node) * np.sin(lat_arg) + \
        np.cos(ascending_node) * np.cos(lat_arg) * np.cos(inclination)
    matrix32 = np.cos(lat_arg) * np.sin(inclination)

    matrix13 = np.sin(ascending_node) * np.sin(inclination)
    matrix23 = - np.cos(ascending_node) * np.sin(inclination)
    matrix33 = np.cos(inclination)

    matrix_ICS2OCS = np.array([[matrix11, matrix12, matrix13],
                               [matrix21, matrix22, matrix23],
                               [matrix31, matrix32, matrix33]])

    pos_vec = np.fromiter(pos.values(), dtype=float)
    target_point = np.fromiter(target_point.values(), dtype=float)
    targ_vec = pos_vec - target_point
    unit_vectors = np.cross(pos_vec, targ_vec) / \
        np.linalg.norm(np.cross(pos_vec, targ_vec), ord=2)
    tetta = np.arccos(np.dot(targ_vec, pos_vec) /
                      (np.linalg.norm(targ_vec, ord=2) * np.linalg.norm(pos_vec, ord=2)))

    q_ICS2OCS = pyquaternion.Quaternion(matrix=matrix_ICS2OCS)
    q_OCS2FCS = pyquaternion.Quaternion(axis=unit_vectors, angle=tetta)
    target_q = q_ICS2OCS * q_OCS2FCS

    return {'w': target_q[0], 'lambda1': target_q[1], 'lambda2': target_q[2], 'lambda3': target_q[3]}


class Satellite():
    """
    класс описывает параметры движения спутника
    """
    def __init__(self, start_pos, start_vel, time_begin, target_point):
        """
        Args:
            start_pos (dict): словарь, содержащий начальное положение спутника
            в ИСК, вида {'x': x, 'y': y,'z': z}, где
                x (float): x - компонента радиус-вектора спутника 
                y (float): y - компонента радиус-вектора спутника 
                z (float): z - компонента радиус-вектора спутника 
            start_vel (dict): словарь, содержащий начальный вектор \
            скорости спутника в ИСК, вида {'vel_x': vel_x, 'vel_y': vel_y,
            'vel_z': vel_z}, где 
                vel_x (float): x - компонента скорости спутника 
                vel_y (float): y - компонента скорости спутника 
                vel_z (float): z - компонента скорости спутника 
            time_begin (float): начальный момент времени
            target_point (dict) : координаты целевой точки наведения в ИСК, 
            вида {'x':x, 'y':y, 'z': z}, где
                x (float): x - координата целевой точки
                y (float): y - координата целевой точки
                z (float): z - координата целевой точки
        """
        self._orb_params = construct_orbital_params(start_pos, start_vel)
        self._transf_mat = None
        self._time_begin = time_begin
        self._current_time = time_begin
        self.pos = {'x': None, 'y': None, 'z': None}
        self.vel = {'vel_x': None, 'vel_y': None, 'vel_z': None}
        self.quat = {'w': None, 'lambda1': None, 'lambda2': None}
        self._target_point = target_point

    def update_params(self, time):
        """
        расчет векторов положения и скорости спутника на момент времени time
        Args:
            time (float): время в секундах  
        """
        self._orb_params.update_anomaly(time)
        self._current_time += time
        self.pos = count_position(self._orb_params._semi_major_axis,
            self._orb_params._ecc_anomaly, self._orb_params._eccentricity, self._transf_mat)

        self.vel = count_velocity(self._orb_params._semi_major_axis,
            self._orb_params._ecc_anomaly, self._orb_params._eccentricity, self._transf_mat)

        self.quat = count_quaternion(self._orb_params._true_anomaly,
            self._orb_params._periapsis_arg, self._orb_params._inclination,
            self._orb_params._ascending_node, self.pos, self._target_point)

    def get_params(self):
        """
        возврат вектора состояния спутника
        Returns:
            result (list): вектор состояния спутника в момент времени time вида
            [julian_date, x, y, z, vel_x, vel_y, vel_z], где
            julian_date (float): время в формате количества юлинский дней в шкале TT, прошедших от эпохи J2000 
            x (float): x - компонента радиус-вектора спутника в осях ИСК
            y (float): y - компонента радиус-вектора спутника в осях ИСК
            z (float): z - компонента радиус-вектора спутника в осях ИСК
            vel_x (float): x - компонента скорости спутника в осях ИСК
            vel_y (float): y - компонента скорости спутника в осях ИСК
            vel_z (float): z - компонента скорости спутника в осях ИСК
            w (float): скалярная часть квантерниона 
            lambda1 (float): x - компонента векторной части кватерниона
            lambda2 (float): y - компонента векторной части кватерниона
            lambda3 (float): z - компонента векторной части кватерниона
        """
        julian_date = JD_START + \
            (self._time_begin + self._current_time) / SECONDS_PER_DAY
        result = [julian_date, self.pos['x'], self.pos['y'], self.pos['z'],
            self.vel['vel_x'], self.vel['vel_y'], self.vel['vel_z'], self.quat['w'],
            self.quat['lambda1'], self.quat['lambda2'], self.quat['lambda3']]

        return result


def construct_satellite(start_pos, start_vel, time_begin, target_point):
    """
    создание и инициализация объекта Satellite
    Args:
       start_pos (dict): словарь, содержащий начальное положение спутника
            в ИСК, вида {'x': x, 'y': y,'z': z}, где
                x (float): x - компонента радиус-вектора спутника 
                y (float): y - компонента радиус-вектора спутника 
                z (float): z - компонента радиус-вектора спутника 
        start_velocity (list): словарь, содержащий начальный вектор скорости 
        спутника в ИСК, вида {'vel_x': vel_x, 'vel_y': vel_y,'vel_z': vel_z}, 
        где 
            vel_x (float): x - компонента скорости спутника
            vel_y (float): y - компонента скорости спутника
            vel_z (float): z - компонента скорости спутника
        time_begin (float): начальный момент времени
    Returns:
        satellite(obj): экземляр класса Satellite
    """
    satellite = Satellite(start_pos, start_vel, time_begin, target_point)
    ascending_node = satellite._orb_params._ascending_node
    periapsis_arg = satellite._orb_params._periapsis_arg
    inclination = satellite._orb_params._inclination

    m11 = np.cos(ascending_node) * np.cos(periapsis_arg) - \
        np.sin(ascending_node) * np.sin(periapsis_arg) * np.cos(inclination)
    m21 = np.sin(ascending_node) * np.cos(periapsis_arg) + \
        np.cos(ascending_node) * np.sin(periapsis_arg) * np.cos(inclination)
    m31 = np.sin(periapsis_arg) * np.sin(inclination)

    m12 = - np.cos(ascending_node) * np.sin(periapsis_arg) - \
        np.sin(ascending_node) * np.cos(periapsis_arg) * np.cos(inclination)
    m22 = - np.sin(ascending_node) * np.sin(periapsis_arg) + \
        np.cos(ascending_node) * np.cos(periapsis_arg) * np.cos(inclination)
    m32 = np.cos(periapsis_arg) * np.sin(inclination)

    satellite._transf_mat = np.array([[m11, m21, m31], [m12, m22, m32]])

    return satellite
