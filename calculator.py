import time

import numpy as np
import pandas as pd

from numba import njit
from io import StringIO
from bokeh.layouts import column, layout
from scipy.constants import speed_of_light
from cartography import CoordLatLon, CoordUtm, CoordEcef
from bokeh.tile_providers import ESRI_IMAGERY, get_provider
from bokeh.plotting import figure, curdoc, show, gridplot, output_file
from bokeh.models import LinearColorMapper, ColorBar, ColumnDataSource, Button
from bokeh.models.widgets import DataTable, DateFormatter, TableColumn, TextAreaInput


@njit
def _sphdot(A, B, C, theta, phi):
    i = A * np.sin(phi) * np.cos(theta)
    j = B * np.sin(phi) * np.sin(theta)
    k = C * np.cos(phi)

    return i + j + k


@njit
def _normal_dist(x , mean , sd):
    prob_density = np.exp(-0.5*((x-mean)/sd)**2)
    return prob_density


@njit
def _beam(theta, phi, theta_d, phi_d, lambda_m):
    fwhm = 1.22 * lambda_m / 6.1
    sig = fwhm / 2.355
    att_phi = _normal_dist(phi, phi_d, sig)
    att_theta = _normal_dist(theta, theta_d, sig)
    return att_phi * att_theta


@njit
def _accumulate(nX, nY, nZ, nants, theta_d, theta, phi_d, phi, lambda_m):
    beam = _beam(theta, phi, theta_d, phi_d, lambda_m)

    tmp = 0 + 0j
    for iant in range(nants):
        tmp += np.exp(
            -1j * 2.0 * np.pi * (
                _sphdot(nX[iant], nY[iant], nZ[iant], theta_d, phi_d) -
                _sphdot(nX[iant], nY[iant], nZ[iant], theta, phi)
            )
        )
    return tmp * beam


class Beams:
    _phi_d: float
    _theta_d: float
    _frequency: float
    _coordinates: pd.DataFrame
    _commited: bool = False

    def __init__(self, frequency: float, phi_d: float,
            theta_d: float, coordinates: pd.DataFrame):
        self.phi_d = phi_d
        self.theta_d = theta_d
        self.frequency = frequency
        self.coordinates = coordinates
        self._commited = True
        self._calculate()

    @property
    def phi_d(self):
        return np.rad2deg(self._phi_d)

    @phi_d.setter
    def phi_d(self, phi_d: float):
        phi_d = float(phi_d)

        if phi_d > 90 or phi_d < 0:
            raise ValueError("Error, zenith angle should be between 0 and 90.")

        self._phi_d = np.deg2rad(phi_d)

        if self._commited:
            self._calculate()

    @property
    def theta_d(self):
        return np.rad2deg(self._theta_d)

    @theta_d.setter
    def theta_d(self, theta_d: float):
        theta_d = float(theta_d)

        if theta_d > 360 or theta_d < 0:
            raise ValueError("Error, azimuth angle should be between 0 and 360.")

        self._theta_d = np.deg2rad(theta_d)

        if self._commited:
            self._calculate()

    @property
    def frequency(self):
        return self._frequency

    @frequency.setter
    def frequency(self, freq: float):
        freq = float(freq) * 1e6

        if freq > 25000e6 or freq < 0:
            raise ValueError("Error, frequency should be between 0 and 25 GHz.")

        self._frequency = freq

        if self._commited:
            self._calculate()

    @property
    def coordinates(self):
        return self._coordinates

    @coordinates.setter
    def coordinates(self, coordinates: pd.DataFrame):
        self._coordinates = coordinates

    @property
    def zenith(self):
        return self._phis

    @property
    def azimuth(self):
        return self._thetas

    @property
    def db(self):
        return self._afdb

    def _calculate(self):
        lambda_m = speed_of_light / self._frequency
        fwhm = 1.22 * lambda_m / 6.1
        sig = fwhm / 2.355

        self._phis = np.linspace(self._phi_d - (sig / 2),
                self._phi_d + (sig / 2), 1024)
        self._thetas = np.linspace(self._theta_d - (sig / 2),
                self._theta_d + (sig / 2), 1024)

        AF = np.zeros(shape=(len(self._thetas), len(self._phis)), dtype=np.complex128)

        nX = np.array(self._coordinates.e / lambda_m, dtype=np.float64)
        nY = np.array(self._coordinates.n / lambda_m, dtype=np.float64)
        nZ = np.array(self._coordinates.u / lambda_m, dtype=np.float64)

        start = time.time()
        for itheta,theta in enumerate(self._thetas):
            for iphi,phi in enumerate(self._phis):
                AF[itheta, iphi] = _accumulate(nX, nY, nZ,
                        len(nX), self._theta_d, theta, self._phi_d, phi, lambda_m)
        print(f"Calculation took {time.time() - start}")

        self._afdb = 10 * np.log10((AF.real * AF.real) + (AF.imag * AF.imag))
        self._afdb -= np.max(self._afdb)


#
# Feed Map
#

antennas_plot = figure(
    title="Array Antennas",
    x_axis_label="East [m]",
    y_axis_label="North [m]",
    x_axis_type="mercator",
    y_axis_type="mercator"
)

color_mapper = LinearColorMapper(palette='Viridis256')

plot = antennas_plot.circle(
    x="lat",
    y="lon",
    size=10.0,
    color={"field": "alt", "transform": color_mapper}
)
antenna_plot_src = plot.data_source

antennas_plot.add_layout(
    ColorBar(
        color_mapper=color_mapper,
        title="Up [m]"
    ),
    place='below'
)

tile_provider = get_provider(ESRI_IMAGERY)
antennas_plot.add_tile(tile_provider)

"""
beam_plot = figure(
    title="Beam Pattern",
    x_axis_label="Azimuth [deg]",
    y_axis_label="Zenith Angle [deg]")

l, r, b, t = np.rad2deg([
    beams.azimuth[-1],
    beams.azimuth[0],
    beams.zenith[-1],
    beams.zenith[0]
])

color_mapper = LinearColorMapper(palette='Viridis256', low=-30.0)

beam_plot.image(
    image=[np.flipud(beams.db.T)],
    x=r, y=b, dw=l-r, dh=b-t,
    color_mapper=color_mapper)

beam_plot.add_layout(
    ColorBar(
        color_mapper=color_mapper,
        title="Beam Gain [dB]"),
    place='below')

azimuth_plot = figure(
    title=f"Azimuth Offset From Az={np.rad2deg(theta_d):.2f}°",
    x_axis_label="Azimuth Offset [deg]",
    y_axis_label="Power [dB]")

azimuth_plot.line(
    x=np.rad2deg(beams.azimuth - beams.theta_d),
    y=beams.db[:, beams.db.shape[1]//2])

zenith_plot = figure(
    title=f"Zenith Angle Offset From ZA={np.rad2deg(phi_d):.2f}°",
    x_axis_label="Zenith Angle Offset [deg]",
    y_axis_label="Power [dB]")

zenith_plot.line(
    x=np.rad2deg(beams.zenith - beams.phi_d),
    y=beams.db[beams.db.shape[0]//2])
"""

with open('ata-datum-condensed.csv', 'r') as file:
    default_ecef = file.read()

ecef_text = TextAreaInput(value=default_ecef, rows=20, title="ECEF CSV File")
ecef_table = DataTable(columns=[
    TableColumn(title="Antenna", field="antenna"),
    TableColumn(title="X", field="x"),
    TableColumn(title="Y", field="y"),
    TableColumn(title="Z", field="z")
])
ecef_table_src = ecef_table.source

def generate_data():
    freq = 500
    phi_d = 30
    theta_d = 30

    ecef_coords = pd.read_csv(StringIO(ecef_text.value))
    latlon_coords = pd.DataFrame(columns=['antenna', 'lat', 'lon', 'alt'])
    enu_coords = pd.DataFrame(columns=['antenna', 'e', 'n', 'u'])

    ref_ant = ecef_coords.iloc[0]
    ref_utm = CoordEcef(ref_ant.x, ref_ant.y, ref_ant.z).to_utm()

    for index, row in ecef_coords.iterrows():
        ant_utm = CoordEcef(row.x, row.y, row.z).to_utm()

        pos = ant_utm.to_latlon().to_webmercator()
        coord = pd.DataFrame({
            "antenna": row.antenna,
            "lat": pos.lat,
            "lon": pos.lon,
            "alt": pos.alt,
        }, index=[0])
        latlon_coords = pd.concat([latlon_coords, coord])

        pos = ant_utm - ref_utm
        coord = pd.DataFrame({
            "antenna": row.antenna,
            "e": pos.e,
            "n": pos.n,
            "u": pos.u,
        }, index=[0])
        enu_coords = pd.concat([enu_coords, coord])

    beams = Beams(freq, phi_d, theta_d, enu_coords)

    ecef_table_src.data = dict(ColumnDataSource(data=ecef_coords).data)
    antenna_plot_src.data = dict(ColumnDataSource(data=latlon_coords).data)


button = Button(label="Press Me")
button.on_click(generate_data)

generate_data()

plots = layout([
    [ecef_text, ecef_table],
    [button],
    [antennas_plot],
], sizing_mode="scale_width")

curdoc().add_root(plots)
