from satellite import construct_satellite
from satellite import JD_START
from main import START_POS, START_VEL
import math 

def test_answer():
    params_time_begin = [JD_START, START_POS['x'], START_POS['y'], START_POS['z'], START_VEL['vel_x'], START_VEL['vel_y'], START_VEL['vel_z']]
    satellite = construct_satellite(START_POS, START_VEL, 0)
    
    params =  satellite.get_params(0)
    compare_flag = True
    for i in range(len(params_time_begin)):
        if not math.isclose(params_time_begin[i], params[i], rel_tol=1e-4):
            compare_flag = False
    
    assert compare_flag
