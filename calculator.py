import utm
import time

import numpy as np
import pandas as pd

from numba import njit
from pyproj import Transformer

from bokeh.plotting import figure, show, gridplot, output_file
from bokeh.models import LinearColorMapper, ColorBar, ColumnDataSource
from bokeh.tile_providers import ESRI_IMAGERY, get_provider

output_file("docs/index.html")

TRAN_4326_TO_3857 = Transformer.from_crs("EPSG:4326", "EPSG:3857")

@njit
def sph2cart(rho, theta, phi):
    x = rho * np.cos(theta) * np.sin(phi)
    y = rho * np.sin(theta) * np.sin(phi)
    z = rho * np.cos(phi)

    return x,y,z

@njit
def sphdot(A, B, C, theta, phi):
    i = A * np.sin(phi) * np.cos(theta)
    j = B * np.sin(phi) * np.sin(theta)
    k = C * np.cos(phi)

    return i + j + k


@njit
def normal_dist(x , mean , sd):
    prob_density = np.exp(-0.5*((x-mean)/sd)**2)
    return prob_density
    #return 1


@njit
def ATABeam(theta, phi, theta_d, phi_d, lambda_m):
    fwhm = 1.22 * lambda_m / 6.1
    sig = fwhm / 2.355
    att_phi = normal_dist(phi, phi_d, sig)
    att_theta = normal_dist(theta, theta_d, sig)
    return att_phi * att_theta


freq = 500#input("Please input observing frequency (in MHz): ") #6500e6 #Hz
freq = float(freq)*1e6
if freq > 11000e6 or freq < 500e6:
    print("Error, observing frequency should be between 500 and 11000 MHz")
    freq = 1400e6

# Define the phase steering angle
#phi_d   = np.deg2rad(30) #np.deg2rad(30)  # Zenith angle
#theta_d = np.deg2rad(0) #np.pi/4  # Azimuth
phi_d = 30#input("Please input the zenith angle [deg]: ")
theta_d = 30#input("Please input the azimuth angle [deg]: ")
phi_d = float(phi_d)
theta_d = float(theta_d)
if phi_d > 90 or phi_d < 0:
    print("Error, zenith angle should be between 0 and 90")
    phi_d = 30

if theta_d > 360 or theta_d < 0:
    print("Error, azimuth should be between 0 and 360")
    theta_d = 30

phi_d = np.deg2rad(phi_d)
theta_d = np.deg2rad(theta_d)



c = 3*10**8 #m/s
lambda_m = c / freq
fwhm = 1.22 * lambda_m / 6.1
sig = fwhm / 2.355

data = pd.read_csv("https://github.com/luigifcruz/ata-datum/raw/main/output_ata_datum.csv")
nX = data["e"] / lambda_m
nY = data["n"] / lambda_m
nZ = data["u"] / lambda_m

nX = np.array(nX, dtype=np.float64)
nY = np.array(nY, dtype=np.float64)
nZ = np.array(nZ, dtype=np.float64)


antennas_plot = figure(
    title="Array Antennas",
    x_axis_label="East [m]",
    y_axis_label="North [m]",
    x_axis_type="mercator",
    y_axis_type="mercator")

ref_data = data.loc[data['antenna'] == "2a"]
ref_lat = float(ref_data["lat_wgs"])
ref_lon = float(ref_data["lon_wgs"])
ref_alt = float(ref_data["alt_wgs"])
ref_x, ref_y, a, b = utm.from_latlon(ref_lat, ref_lon)

web_mercator = {
    "X": [],
    "Y": [],
    "Z": [],
}

for index, row in data.iterrows():
    e = row["e"] + ref_x
    n = row["n"] + ref_y
    u = row["u"] + ref_alt
    ref_lat, ref_lon = utm.to_latlon(e, n, a, b)

    x, y = TRAN_4326_TO_3857.transform(ref_lat, ref_lon)

    web_mercator["X"].append(x)
    web_mercator["Y"].append(y)
    web_mercator["Z"].append(u)

source = ColumnDataSource(data=web_mercator)

color_mapper = LinearColorMapper(palette='Viridis256')

antennas_plot.circle(
    x="X",
    y="Y",
    size=10.0,
    source=source,
    color={"field": "Z", "transform": color_mapper})

antennas_plot.add_layout(
    ColorBar(
        color_mapper=color_mapper,
        title="Up [m]"),
    place='below')

tile_provider = get_provider(ESRI_IMAGERY)
antennas_plot.add_tile(tile_provider)

assert len(nX) == len(nY) == len(nZ)
Nants = len(nX)
print(Nants)

phis   = np.linspace(phi_d - sig/2, phi_d + sig/2, 1024)        # Zenith angle
thetas = np.linspace(theta_d - sig/2, theta_d + sig/2, 1024)    # Azimuth

AF = np.zeros(shape=(len(thetas), len(phis)), dtype=np.complex128)

@njit
def accumulate(nX, nY, nZ, theta_d, theta, phi_d, phi, lambda_m):
    beam = ATABeam(theta, phi, theta_d, phi_d, lambda_m)

    tmp = 0 + 0j
    for iant in range(Nants):
        tmp += np.exp(-1j*2*np.pi * (sphdot(nX[iant], nY[iant], nZ[iant], theta_d, phi_d) -
                                     sphdot(nX[iant], nY[iant], nZ[iant], theta, phi)))
    return tmp * beam


start = time.time()
for itheta,theta in enumerate(thetas):
    for iphi,phi in enumerate(phis):
        AF[itheta, iphi] = accumulate(nX, nY, nZ, theta_d, theta, phi_d, phi, lambda_m)
print(f"Calculation took {time.time() - start}")

AF_db = 10*np.log10(AF.real*AF.real + AF.imag*AF.imag)
print(np.max(AF_db))
AF_db -= np.max(AF_db)

beam_plot = figure(
    title="Beam Pattern",
    x_axis_label="Azimuth [deg]",
    y_axis_label="Zenith Angle [deg]")

l, r, b, t = np.rad2deg([thetas[-1], thetas[0], phis[-1], phis[0]])

color_mapper = LinearColorMapper(palette='Viridis256', low=-30.0)

beam_plot.image(
    image=[np.flipud(AF_db.T)],
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
    x=np.rad2deg(thetas - theta_d),
    y=AF_db[:, AF_db.shape[1]//2])

zenith_plot = figure(
    title=f"Zenith Angle Offset From ZA={np.rad2deg(phi_d):.2f}°",
    x_axis_label="Zenith Angle Offset [deg]",
    y_axis_label="Power [dB]")

zenith_plot.line(
    x=np.rad2deg(phis - phi_d),
    y=AF_db[AF_db.shape[0]//2])

plots = gridplot([
    [antennas_plot, beam_plot],
    [azimuth_plot, zenith_plot],
])
show(plots)
