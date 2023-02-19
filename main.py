from satellite import construct_satellite
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# начальный момент времени
TIME_BEGIN = 0 
# конечный момент времени
TIME_END = 20000
# временной шаг
STEP = 1
# начальное положение спутника в ИСК
START_POS = {'x': 4362521.19692133, 'y': -
             2174459.71448059, 'z': 4720847.40402189}
# начальный вектор скорости спутника в ИСК          
START_VEL = {'vel_x': 5356.39915538069,
             'vel_y': 4741.41348686709, 'vel_z': -2761.5472632395}

# coordinates_X = []
# coordinates_Y = []
# coordinates_Z = []
# velocity_X = []
# velocity_Y = []
# velocity_Z = []
#Ecc = []

def process(time_begin, time_end, step=1):
    """
    запись в файл вектора состояния спутника в заданном промежутке времени
    Args:
        time_begin (float): начальный момент времени
        time_end (_type_): конечный момент времени
        step (int, optional): временной шаг
    """
    satellite = construct_satellite(START_POS, START_VEL, TIME_BEGIN)
    with open('result.txt', 'w') as f:
        for i in range(time_begin, time_end, step):
            out_data = satellite.get_params(i)
            
            # Ecc.append(satellite._orb_params._ecc_anomaly)
            # coordinates_X.append(out_data[1])
            # coordinates_Y.append(out_data[2])
            # coordinates_Z.append(out_data[3])
            # velocity_X.append(float(vel_x))
            # velocity_Y.append(float(vel_y))
            # velocity_Z.append(float(vel_z))  
            out_data = ' '.join(str(i) for i in out_data) + '\n'
            f.writelines(out_data)

if __name__ == "__main__":
    process(TIME_BEGIN, TIME_END, STEP)

#     # print(E)
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
#     ax.plot(coordinates_X, coordinates_Y, coordinates_Z)
#     #ax.plot(vel_x, coordinates_Y, coordinates_Z)
#    # plt.plot((coordinates_X, coordinates_Y)
# plt.plot(Ecc)
# plt.show()
