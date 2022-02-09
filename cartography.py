import utm
from pyproj import Transformer


_latlon_wgs_proj = {"proj": 'latlong', "ellps": 'WGS84', "datum": 'WGS84'}
_ecef_wgs_proj = {"proj": 'geocent', "ellps": 'WGS84', "datum": 'WGS84'}
_latlon_webmercator_proj = "EPSG:3857"


class CoordLatLon:
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

    def __str__(self) -> str:
        return f"""Lat/Lon WGS84 Coordinate:
    Latitude  (deg): {self.lat:.6f}
    Longitude (deg): {self.lon:.6f}
    Altitude    (m): {self.alt:.6f}"""

    def to_webmercator(self):
        recipe = Transformer.from_crs(_latlon_wgs_proj, _latlon_webmercator_proj)
        lat, lon, alt = recipe.transform(self.lon, self.lat, self.alt)
        return CoordLatLon(lat, lon, alt)

    def to_ecef(self):
        recipe = Transformer.from_crs(_latlon_wgs_proj, _ecef_wgs_proj)
        x, y, z = recipe.transform(self.lon, self.lat, self.alt)
        return CoordEcef(x, y, z)

    def to_utm(self):
        e, n, zone_number, zone_letter = utm.from_latlon(self.lat, self.lon)
        return CoordUtm(e, n, self._alt, zone_number, zone_letter)


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

    def __str__(self) -> str:
        return f"""ECEF Coordinate:
    X (m): {self.x:.6f}
    Y (m): {self.y:.6f}
    Z (m): {self.z:.6f}"""

    def to_latlon(self):
        recipe = Transformer.from_crs(_ecef_wgs_proj, _latlon_wgs_proj)
        lon, lat, alt = recipe.transform(self.x, self.y, self.z)
        return CoordLatLon(lat, lon, alt)

    def to_utm(self):
        return self.to_latlon().to_utm()


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
        return self._e

    @property
    def n(self) -> float:
        return self._n

    @property
    def u(self) -> float:
        return self._u

    @property
    def zone_number(self) -> int:
        return self._zone_number

    @property
    def zone_letter(self) -> str:
        return self._zone_letter

    def __str__(self) -> str:
        return f"""UTM Coordinate:
    Easting  (m): {self.e:.6f}
    Northing (m): {self.n:.6f}
    Up       (m): {self.u:.6f}
    Zone        : {self.zone_number}{self.zone_letter}"""

    def to_latlon(self):
        lat, lon = utm.to_latlon(self.e, self.n, self.zone_number, self.zone_letter)
        return CoordLatLon(lat, lon, self.u)

    def to_ecef(self):
        return self.to_latlon().to_ecef()

    def __add__(self, other):
        if self._zone_number != other._zone_number:
            raise ValueError("Incompatible zone number.")

        if self._zone_letter != other._zone_letter:
            raise ValueError("Incompatible zone letter.")

        e = self.e + other.e
        n = self.n + other.n
        u = self.u + other.u

        return CoordUtm(e, n, u, self.zone_number, self.zone_letter)

    def __sub__(self, other):
        if self._zone_number != other._zone_number:
            raise ValueError("Incompatible zone number.")

        if self._zone_letter != other._zone_letter:
            raise ValueError("Incompatible zone letter.")

        e = self.e - other.e
        n = self.n - other.n
        u = self.u - other.u

        return CoordUtm(e, n, u, self.zone_number, self.zone_letter)
