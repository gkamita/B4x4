#!/usr/bin/python
# encoding: utf-8
"""
B4x4
====

Provides
  An interface for simulating, plotting and exporting
  reflectance/transmittance spectra of chiral nematic stacks
  by the Berreman 4x4 matrix method.
  Berreman4x4.py is required for using this module.

Documentation
  Available as docstrings.

Classes
  Factory: In charge of performing Berreman 4x4 simulation.
           Instantiate Factory to start your simulation.
           When a simulation is performed, Factory returns
           an instance of Spectrum as an output.

  Spectrum: In charge of plotting and exporting the results.

Author
  Gen Kamita (gk298@cam.ac.uk)
"""

__version_ = '1.0.0'

import Berreman4x4
from numpy import sin, abs, array, linspace, append, arange, savetxt
from Berreman4x4 import pi, e_y
from matplotlib import pyplot
from time import localtime
from collections import namedtuple
from copy import deepcopy


class Factory(object):
    """
    Performes simulations using properties stored in class instances.
    The default parameteres are based on Dumanli et. al.ACS Appl. Mater. Interfaces 2014, 6, 12302,
    which is suitable for simulating dry wood pulp CNC films.

    Methods
      calculateL: Returns instance of Spectrum containing LCP reflectance.
      calculateR: Same as calculateL but for RCP reflectance.
      matrix: Run simulations in batches.
    
    Properties
      pitch:        pitch in nm (180 degree twist)
      no:           refractive index that the ordinary ray experiences
      ne:           refractive index that the extraordinary ray experiences
      nAverage:     Average refractive index, that is (no + ne)/2
      nDelta:       Birefringence, that is ne - no
      nSuperstrate: refractive index of superstrate (air: n = 1)
      nSubstrate:   refractive index of suberstrate (PS: n = 1.59, SiO2: n=1.55)
      stack:        number of pitches (180 degree twist) in the structure
      angle:        angle of incidence in degrees
      lbda_min:     wavelength min (in meters)
      lbda_max:     wavelength max (in meters)
      points:       number of data points in simulation
      slices:       number of discrete anisotropic layers per 180 twist

    Attribute
      settings:     All parameters related to the above properties are stored here.
    """

    def __init__(self):
        self.settings = dict(no=1.524, ne=1.586, nSubstrate=1.59, nSuperstrate=1, pitch=190, stack=10, angle=0,
                             lbda_min=400,
                             lbda_max=800, points=101, slices=23)

    def prepare(self):
        sup = Berreman4x4.IsotropicNonDispersiveMaterial(self.settings['nSuperstrate'])
        sub = Berreman4x4.IsotropicNonDispersiveMaterial(self.settings['nSubstrate'])
        front = Berreman4x4.IsotropicHalfSpace(sup)
        back = Berreman4x4.IsotropicHalfSpace(sub)

        LC = Berreman4x4.UniaxialNonDispersiveMaterial(self.settings['no'], self.settings['ne'])  # ne along z
        R = Berreman4x4.rotation_v_theta(e_y, pi / 2)  # rotation of pi/2 along y
        LC = LC.rotated(R)  # apply rotation from z to x
        # Cholesteric pitch:

        # One half turn of a left-handed helix. To change to right-hand, angle=+pi.
        TN = Berreman4x4.TwistedMaterial(LC, self.settings['pitch'], angle=-pi, div=self.settings['slices'])

        # Inhomogeneous layer, repeated layer, and structure
        IL = Berreman4x4.InhomogeneousLayer(TN)
        # h = N * p/2
        L = Berreman4x4.RepeatedLayers([IL], self.settings['stack'])
        s = Berreman4x4.Structure(front, [L], back)
        lbda_list = wavelength(self.settings)
        k0_list = 2 * pi / lbda_list

        # Calculation with Berreman4x4
        J = array([s.getJones(sin(self.settings['angle'] * pi / 180), k0) for k0 in k0_list])
        # power = abs(J)**2

        # Jones matrices for the circular wave basis
        Jc = Berreman4x4.circularJones(J)
        power_c = abs(Jc) ** 2
        return power_c

    def calculate(self):
        power_c = self.prepare()
        # Right-circular wave is reflected in the stop-band.
        # R_LR, T_LR close to zero.
        # r = Spectrum(Berreman4x4.extractCoefficient(power_c, 'r_RR'), self.settings,
        # type='right circularly polarised light')
        # l = Spectrum(Berreman4x4.extractCoefficient(power_c, 'r_LL'), self.settings)
        r = Spectrum(Berreman4x4.extractCoefficient(power_c, 't_pp'), self.settings,
                     type='right circularly polarised light')
        l = Spectrum(Berreman4x4.extractCoefficient(power_c, 't_ss'), self.settings)

        LandR = namedtuple('LandR', 'l r')
        result = LandR(l, r)
        return result

    def calculateR(self):
        power_c = self.prepare()
        # Right-circular wave is reflected in the stop-band.
        # R_LR, T_LR close to zero.
        spectrum = Spectrum(Berreman4x4.extractCoefficient(power_c, 'r_RR'), self.settings, type='r_left')
        return spectrum

    def calculateL(self):
        power_c = self.prepare()
        # Left-circular wave is transmitted in the full spectrum.
        # T_RL, R_RL, R_LL close to zero, T_LL close to 1.
        spectrum = Spectrum(Berreman4x4.extractCoefficient(power_c, 'r_LL'), self.settings)
        # transmission = Spectrum(1 - reflection.spectrum, self.settings, type='t_left')
        # r_and_t = namedtuple('r_and_t', 'r t')
        # spectra = r_and_t(reflection, transmission)
        return spectrum

    def matrix(self, parameter, start, stop, step=1):
        """
        Returns batch-simulated spectra.
        parameter: String that specifis which parameter to vary. 
        """
        leftArray = array([])
        if parameter in self.settings:

            for x in arange(start, stop, step):
                self.settings[parameter] = x
                result = self.calculateL()
                leftArray = append(leftArray, result.spectrum)

        elif parameter is 'nDelta':

            for x in arange(start, stop, step):
                self.nDelta = x
                result = self.calculateL()
                leftArray = append(leftArray, result.spectrum)

        elif parameter is 'nAverage':

            for x in arange(start, stop, step):
                self.nAverage = x
                result = self.calculateL()
                leftArray = append(leftArray, result.spectrum)

        else:
            return 0
        leftMatrix = leftArray.reshape(-1, self.settings['points'])
        Result = Spectrum(leftMatrix, self.settings)
        return Result


    @property
    def nDelta(self):
        return self.settings['ne'] - self.settings['no']

    @nDelta.setter
    def nDelta(self, nDelta):
        half = nDelta / 2.0
        nAverage = self.nAverage
        self.settings['no'] = nAverage - half
        self.settings['ne'] = nAverage + half

    @property
    def nAverage(self):
        return (self.settings['no'] + self.settings['ne']) / 2.0

    @nAverage.setter
    def nAverage(self, nAverage):
        nDelta = self.nDelta
        self.settings['no'] = nAverage - nDelta / 2.0
        self.settings['ne'] = nAverage + nDelta / 2.0

    @property
    def no(self):
        return self.settings['no']

    @no.setter
    def no(self, no):
        self.settings['no'] = no

    @property
    def ne(self):
        return self.settings['ne']

    @ne.setter
    def ne(self, ne):
        self.settings['ne'] = ne

    @property
    def pitch(self):
        return self.settings['pitch']

    @pitch.setter
    def pitch(self, pitch):
        self.settings['pitch'] = pitch

    @property
    def nSuperstrate(self):
        return self.settings['nSuperstrate']

    @nSuperstrate.setter
    def nSuperstrate(self, nSuperstrate):
        self.settings['nSuperstrate'] = nSuperstrate

    @property
    def nSubstrate(self):
        return self.settings['nSubstrate']

    @nSubstrate.setter
    def nSubstrate(self, nSubstrate):
        self.settings['nSubstrate'] = nSubstrate

    @property
    def stack(self):
        return self.settings['stack']

    @stack.setter
    def stack(self, stack):
        self.settings['stack'] = stack

    @property
    def angle(self):
        return self.settings['angle']

    @angle.setter
    def angle(self, angle):
        self.settings['angle'] = angle

    @property
    def lbda_min(self):
        return self.settings['lbda_min']

    @lbda_min.setter
    def lbda_min(self, lbda_min):
        self.settings['lbda_min'] = lbda_min

    @property
    def lbda_max(self):
        return self.settings['lbda_max']

    @lbda_max.setter
    def lbda_max(self, lbda_max):
        self.settings['lbda_max'] = lbda_max

    @property
    def points(self):
        return self.settings['points']

    @points.setter
    def points(self, points):
        self.settings['points'] = points

    @property
    def slices(self):
        return self.settings['slices']

    @slices.setter
    def slices(self, slices):
        self.settings['slices'] = slices


class Spectrum(object):
    """
    Output of Factory.

    Method
      plot: Plots the spectrum.
      save: Exports the spectrum as a txt file in the current directory.
            It will also save the settings and time it was saved in the header.

    Attribute
      spectrum: The simulated spectrum.
      settings: The setting which was used for the simulation are stored here.
    """

    def __init__(self, spectrum, settings, type='r_left'):
        self.spectrum = spectrum
        self.settings = deepcopy(settings)
        self.settings['type'] = type

    def plot(self):
        fig = pyplot.figure()
        ax = fig.add_subplot("111")
        lbda_list = wavelength(self.settings)
        ax.plot(lbda_list, self.spectrum, '--', label=self.settings['type'])
        ax.set_xlabel(r"Wavelength $\lambda_0$ (nm)")
        ax.set_ylabel(r"Reflectance $R$")
        fmt = ax.xaxis.get_major_formatter()
        fmt.set_powerlimits((-3, 3))
        pyplot.show()

    def image(self):
        lbda_list = wavelength(self.settings)
        fig = pyplot.imshow(self.spectrum, extent=[self.settings['lbda_min'], self.settings['lbda_max'], 90, 0])
        pyplot.xlabel(r"Wavelength [nm]")
        pyplot.ylabel(r"Angle [degrees]")
        pyplot.colorbar()
        pyplot.show()


    def save(self, filename='spectrum.txt'):
        self.settings['savedtime'] = localtime()
        try:
            savetxt(filename, self.spectrum, header=str(self.settings))
        except TypeError:
            savetxt(filename, self.spectrum)
            print "Header not included in the output text due to the version of your numpy module."


def wavelength(settings):
    """
    Used for calculating wavelength.
    :param settings: Dictionary "B4x4.settings"
    :return: ndarray of wavelength
    """
    return linspace(settings['lbda_min'], settings['lbda_max'], settings['points'])


if __name__ == '__main__':
    mySim = Factory()
    mySpec = mySim.calculateL()
    mySpec.plot()
    # mySpec.save('mytResult.txt')

