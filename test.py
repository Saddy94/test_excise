from pymap3d import ecef2eci
from pymap3d import geodetic2ecef
import datetime


time = datetime.datetime(2022, 2, 18, 13, 21, 16, 584)
lat = 45.920266
lon =  63.342286

x, y, z = geodetic2ecef(lat, lon, 0)
x_eci, y_eci, z_eci = ecef2eci(x,y,z,time)
print(x,y,z)
print(float(x_eci), float(y_eci), float(z_eci))