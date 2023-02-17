from math import asin, atan, atan2, cos, degrees, pi, radians, sin, sqrt


def llh2xyz(lat: float, lon: float, h, a=6378137.0, finv=298.257223563, input_degrees=True):
    """
    Converts ellipsoidal coordinates to cartesian.
    :param lat:
    :param lon:
    :param h:
    :return: A tuple containing cartesian coordinates in ECEF [X, Y, Z] meters.
    """
    # convert degrees to radians if necessary
    if input_degrees:
        lat = radians(lat)
        lon = radians(lon)
    # Default (a=6378137.0, finv=298.257223563) is for WGS84
    f = 1 / finv
    b = a * (1 - f)
    e2 = 1 - (1 - f) ** 2
    # Compute the Cartesian coordinates
    v = a / sqrt(1 - e2 * sin(lat) * sin(lat))
    x = (v + h) * cos(lat) * cos(lon)
    y = (v + h) * cos(lat) * sin(lon)
    z = (v * (1 - e2) + h) * sin(lat)
    return x, y, z