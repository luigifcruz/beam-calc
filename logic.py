import numpy as np
import pandas as pd

from numba import njit
from scipy.constants import speed_of_light


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
        return np.rad2deg(self._phis)

    @property
    def azimuth(self):
        return np.rad2deg(self._thetas)

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

        for itheta,theta in enumerate(self._thetas):
            for iphi,phi in enumerate(self._phis):
                AF[itheta, iphi] = _accumulate(nX, nY, nZ,
                        len(nX), self._theta_d, theta, self._phi_d, phi, lambda_m)

        self._afdb = 10 * np.log10((np.real(AF)**2) + (np.imag(AF)**2))
        self._afdb -= np.max(self._afdb)
