import utm
from pyproj import Transformer
from functools import wraps


def _add_method(cls):
    def decorator(func):
        @wraps(func)
        def wrapper(self, *args, **kwargs):
            return func(self, *args, **kwargs)
        setattr(cls, func.__name__, wrapper)
        return func
    return decorator


class CoordWgs:
    _lat: float
    _lon: float
    _alt: float

    def __init__(self, lat: float, lon: float, alt: float):
        self._lat = lat
        self._lon = lon
        self._alt = alt

    @property
    def lat(self) -> float:
        return self._lat

    @property
    def lon(self) -> float:
        return self._lon

    @property
    def alt(self) -> float:
        return self._alt


class CoordEcef:
    _x: float
    _y: float
    _z: float

    def __init__(self, x: float, y: float, z: float):
        self._x = x
        self._y = y
        self._z = z

    @property
    def x(self) -> float:
        return self._x

    @property
    def y(self) -> float:
        return self._y

    @property
    def z(self) -> float:
        return self._z


class CoordUtm:
    _e: float
    _n: float
    _u: float
    _zone_number: int
    _zone_letter: str

    def __init__(self, e: float, n: float, u: float,
            zone_number: int, zone_letter):
        self._e = e
        self._n = n
        self._u = u
        self._zone_number = int(zone_number)
        self._zone_letter = zone_letter

    @property
    def e(self) -> float:
        return self._u

    @property
    def n(self) -> float:
        return self._n

    @property
    def u(self) -> float:
        return self._u


_wgs_proj = {"proj": 'latlong', "ellps": 'WGS84', "datum": 'WGS84'}
_ecef_proj = {"proj": 'geocent', "ellps": 'WGS84', "datum": 'WGS84'}


@_add_method(CoordWgs)
def to_ecef(self) -> CoordEcef:
    recipe = Transformer.from_crs(_wgs_proj, _ecef_proj)
    x, y, z = recipe.transform(self._lon, self._lat, self._alt)
    return CoordEcef(x, y, z)


@_add_method(CoordWgs)
def to_utm(self) -> CoordUtm:
    e, n, zone_number, zone_letter = utm.from_latlon(
        self._lat, self._lon, self._zone_number, self._zone_letter)
    return CoordUtm(e, n, self._alt, zone_number, zone_letter)


@_add_method(CoordUtm)
def to_wgs(self) -> CoordWgs:
    lat, lon = utm.to_latlon(self._e, self._n, self._zone_number, self._zone_letter)
    return CoordWgs(lat, lon, self._u)


@_add_method(CoordUtm)
def to_ecef(self) -> CoordEcef:
    return self.to_wgs().to_ecef()


@_add_method(CoordEcef)
def to_wgs(self) -> CoordWgs:
    recipe = Transformer.from_crs(_ecef_proj, _wgs_proj)
    lon, lat, alt = recipe.transform(self._x, self._y, self._z)
    return CoordWgs(lat, lon, alt)


@_add_method(CoordEcef)
def to_utm(self) -> CoordUtm:
    return self.to_wgs().to_utm()
