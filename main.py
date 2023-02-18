from satellite import construct_satellite
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

"""
:param TIME_BEGIN: int, начальный момент времени в секундах
:param TIME_END: int, конечный момент времени в секундах
:param START_POS: list[float], начальное положение ка 
:param START_VEL: list[float], начальные составляющие вектор скорости ка
"""
TIME_BEGIN = 0
TIME_END = 20000
START_POS = {'x': 4362521.19692133, 'y': -
             2174459.71448059, 'z': 4720847.40402189}
START_VEL = {'vel_x': 5356.39915538069,
             'vel_y': 4741.41348686709, 'vel_z': -2761.5472632395}


def convert_data_to_str(data):
    result = ' '.join(str(i) for i in data) + '\n'
    return result

# coordinates_X = []
# coordinates_Y = []
# coordinates_Z = []
# velocity_X = []
# velocity_Y = []
# velocity_Z = []
# Ecc = []

# TODO: разобраться с шагом, чтобы можно было разные шаги делать, а не только delta_t=1
def process(time_begin, time_end, delta_t=1):
    satellite = construct_satellite(START_POS, START_VEL)
    with open('result.txt', 'w') as f:
        for i in range(int(time_end)):
            out_data = satellite.get_params(i)
            # coordinates_X.append(out_data[1])
            # coordinates_Y.append(out_data[2])
            # coordinates_Z.append(out_data[3])
            # velocity_X.append(float(vel_x))
            # velocity_Y.append(float(vel_y))
            # velocity_Z.append(float(vel_z))  
            out_data = convert_data_to_str(out_data)
            f.writelines(out_data)

if __name__ == "__main__":
    process(TIME_BEGIN, TIME_END)

#     # print(E)
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
#     ax.plot(coordinates_X, coordinates_Y, coordinates_Z)
#     #ax.plot(vel_x, coordinates_Y, coordinates_Z)
#    # plt.plot((coordinates_X, coordinates_Y)
#    # plt.plot(Ecc)
#     plt.show()
