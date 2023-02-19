from orbital_params import construct_orbital_params
import numpy as np
import datetime
from pymap3d import ecef2eci, geodetic2ecef
import pyquaternion
from physical_const import *


# начальный момент времени в юлианских днях
JD_START = 8084.05644194318 
# начальный момент времени в формате utc
UTC_TIME_START = datetime.datetime(2022, 2, 18, 13, 21, 16, 584)
# количество секунд в сутках
SECONDS_PER_DAY = 86400
# координаты целевой точки на поверхности Земли
TARGET_POINT = {'lat': 45.920266, 'lon': 63.342286}


class Satellite():
    """
    класс описывает параметры движения спутника
    """
    def __init__(self, start_pos, start_vel, time_begin):
        """
        Args:
            start_pos (dict): словарь, содержащий начальное положение спутника,
            вида {'x': x, 'y': y,'z': z}, где
                x (float): x - компонента радиус-вектора спутника в осях ИСК
                y (float): y - компонента радиус-вектора спутника в осях ИСК
                z (float): z - компонента радиус-вектора спутника в осях ИСК
            start_vel (dict): словарь, содержащий начальный вектор \
            скорости спутника в ИСК, вида {'vel_x': vel_x, 'vel_y': vel_y,
            'vel_z': vel_z}, где 
                vel_x (float): x - компонента скорости спутника в осях ИСК
                vel_y (float): y - компонента скорости спутника в осях ИСК
                vel_z (float): z - компонента скорости спутника в осях ИСК
            time_begin (float): начальный момент времени
        
        _orb_params (obj): кеплеровы параметры орбиты
        _transform_matrix (ndarray) матрица пересчета в ИСК
        _transform_matrix1 (ndarray) матрица пересчета в ИСК
        """
        self._orb_params = construct_orbital_params(start_pos, start_vel)
        self._transform_matrix = None
        self._transform_matrix1 = None
        self._time_begin = time_begin
    
    def _count_position(self):
        """
        расчет положения спутника на текущий момент времени
        Returns:
            x (float): x - компонента радиус-вектора спутника в осях ИСК
            y (float): y - компонента радиус-вектора спутника в осях ИСК
            z (float): z - компонента радиус-вектора спутника в осях ИСК
        """
        x, y, z = (self._orb_params._semi_major_axis * (self._transform_matrix * (np.cos(self._orb_params._ecc_anomaly) - self._orb_params._eccentricity) + self._transform_matrix1 *
                                                        np.sqrt(1 - self._orb_params._eccentricity ** 2) * np.sin(self._orb_params._ecc_anomaly)))
        return float(x), float(y), float(z)

    def _count_velocity(self):
        """
        расчет вектора скорости спутника на текущий момент времени
        Returns:
            vel_x (float): x - компонента скорости спутника в осях ИСК
            vel_y (float): y - компонента скорости спутника в осях ИСК
            vel_z (float): z - компонента скорости спутника в осях ИСК
        """
        vel_x, vel_y, vel_z = (np.sqrt(MU / self._orb_params._semi_major_axis) / (1 - self._orb_params._eccentricity * np.cos(self._orb_params._ecc_anomaly))) * ((- self._transform_matrix) * np.sin(self._orb_params._ecc_anomaly) + self._transform_matrix1 * (
            np.sqrt(1 - self._orb_params._eccentricity ** 2) * np.cos(self._orb_params._ecc_anomaly)))

        return float(vel_x), float(vel_y), float(vel_z)

    def count_quanternion(self):
        """
        расчет квантерниона целевой ориентации на текущий момент времени
        """
        pass

    def get_params(self, time):
        """
        расчет векторов положения и скорости спутника на момент времени time
        Args:
            time (float): время в секундах

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
        self._orb_params.update_anomaly(time)

        julian_date = JD_START + (self._time_begin + time) / SECONDS_PER_DAY 

        x, y, z = self._count_position()
        vel_x, vel_y, vel_z = self._count_velocity()
        result = [julian_date, x, y, z, vel_x, vel_y, vel_z]
        
        return result


def construct_satellite(start_pos, start_vel, time_begin):
    """
    создание и инициализация объекта Satellite
    Args:
        start_pos (dict): словарь, содержащий начальное положение спутника,
        вида {'x': x, 'y': y,'z': z}, где
            x (float): x - компонента радиус-вектора спутника в осях ИСК
            y (float): y - компонента радиус-вектора спутника в осях ИСК
            z (float): z - компонента радиус-вектора спутника в осях ИСК
        start_velocity (list): словарь, содержащий начальный вектор скорости 
        спутника в ИСК, вида {'vel_x': vel_x, 'vel_y': vel_y,'vel_z': vel_z}, 
        где 
            vel_x (float): x - компонента скорости спутника в осях ИСК
            vel_y (float): y - компонента скорости спутника в осях ИСК
            vel_z (float): z - компонента скорости спутника в осях ИСК
        time_begin (float): начальный момент времени
    Returns:
        satellite(obj): экземляр класса Satellite
    """
    satellite = Satellite(start_pos, start_vel, time_begin)

    ascending_node = satellite._orb_params._ascending_node
    periapsis_arg = satellite._orb_params._periapsis_arg
    inclination = satellite._orb_params._inclination
    #TODO переделать формат матрицы
    satellite._transform_matrix = np.array([[np.cos(ascending_node) * np.cos(periapsis_arg) - np.sin(ascending_node) * np.sin(periapsis_arg) * np.cos(inclination)],
                                            [np.sin(ascending_node) * np.cos(periapsis_arg) + np.cos(ascending_node)
                                             * np.sin(periapsis_arg) * np.cos(inclination)],
                                            [np.sin(periapsis_arg) * np.sin(inclination)]])
    satellite._transform_matrix1 = np.array([[- np.cos(ascending_node) * np.sin(periapsis_arg) - np.sin(ascending_node) * np.cos(periapsis_arg) * np.cos(inclination)],
                                             [-np.sin(ascending_node) * np.sin(periapsis_arg) + np.cos(ascending_node)
                                              * np.cos(periapsis_arg) * np.cos(inclination)],
                                             [np.cos(periapsis_arg) * np.sin(inclination)]])

    return satellite


def convert_geo2eci(target_point, time):
    """
    преобразование координат целевой точки из геодезической СК в ИСК
    Args:
        target_point (dict): координатцы целевой точки в геодезической СК
        time (datetime): время наблюдения

    Returns:
        eci_coord (list): список, содержащий координаты целевой точки в ИСК вида
        [x, y, z], где
            x (float): x - компонента радиус-вектора целевой точки в осях ИСК
            y (float): y - компонента радиус-вектора целевой точки в осях ИСК
            z (float): z - компонента радиус-вектора целевой точки в осях ИСК
    """
    x, y, z = geodetic2ecef(target_point['lat'], target_point['lon'], 0)
    x, y, z = ecef2eci(x, y, z, time)
    eci_coord = [x, y, z]
    return eci_coord
