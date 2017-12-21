#!/usr/bin/env python
"""
Compute the statistics of sampled sets produced by OpenFOAM
"""

from functools import wraps
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quad
import warnings

from HbO2.setup.simulationParameters import IOHbO2ParametersAxisymmetric
from postprocessing.loadSampledSetFiles import loadSampledSets


class SampledSetStatError(Exception): pass
class TooLowAmplitudeError(SampledSetStatError): pass
class TooHighAmplitudeError(SampledSetStatError): pass


def catchAmplitudeError(lowAmplValue, highAmplValue):
    """Catch amplitude errors and sets given values if errors are caught.

    Args:
        lowAmplValue:  return value set if TooLowAmplitudeError is caught
        highAmplValue: return value set if TooHighAmplitudeError is caught

    Returns:
        Wrapper function
    """
    def wrap(f):
        @wraps(f)
        def wrapped_f(*args):
            try:
                r = f(*args)
            except TooLowAmplitudeError:
                warnings.warn('Too low amplitudes, setting amplitude radius to %g.' 
                                 % lowAmplValue)
                r = lowAmplValue
            except TooHighAmplitudeError:
                warnings.warn('Too high amplitudes, setting amplitude radius to %g.'
                                 % highAmplValue)
                r = highAmplValue
            return r
        return wrapped_f
    return wrap


class SampledSetStatistics:

    def __init__(self, sampledSetsDict):
        self.sampledSetsDict = sampledSetsDict

    def positions(self):
        return self.sampledSetsDict['positions']

    def average(self, minTime=0.0):
        """Compute the time-averaged value of the sampled set."""
        timeKeys = self.timeKeys(minTime)
        return np.average(np.asarray([self.sampledSetsDict[t] for t in timeKeys]),
                          axis=0)

    def fluctuationAmplitudes(self, timeWindow=0.0):
        """Compute the fluctuation amplitude of the sampled set at the latest times.

        Args: 
            timeWindow: time over which the amplitude is computed
        Returns:
            1D numpy array with amplitudes at the loaded positions
        """
        minima, maxima = self.lastMinMax(timeWindow)
        return maxima - minima

    def interpolatedAmplitude(self, r, timeWindow=0.0):
        """Compute the fluctuation amplitude at a given radius.

        Args:
            r: radius where to compute the amplitude
            timeWindow: time over which the amplitude is computed
        Returns:
            fluctuation amplitude
        """
        f = interp1d(self.positions(), self.fluctuationAmplitudes(timeWindow), kind='cubic')
        return f(r)

    def interpolatedRelativeAmplitude(self, r, timeWindow=0.0):
        """Compute the relative fluctuation amplitude at a given radius.

        Args:
            r: radius where to compute the amplitude
            timeWindow: time over which the amplitude is computed
        Returns:
            relative fluctuation amplitude
        """
        amplitude = self.interpolatedAmplitude(r, timeWindow)
        maxAmplitude = np.max(self.fluctuationAmplitudes(timeWindow))
        return amplitude/maxAmplitude

    def amplitudePosition(self, amplitude, timeWindow):
        """Compute the position at which the given amplitude is reached.
        """
        x = self.positions()
        a = self.fluctuationAmplitudes(timeWindow)
        try:
            idx = np.where(a >= amplitude)[0][-1]
        except IndexError:
            raise TooLowAmplitudeError, \
                "No element with amplitude greater than %g was found." % amplitude
        
        try: 
            return x[idx] + (amplitude - a[idx]) \
                          * (x[idx+1] - x[idx])/(a[idx+1] - a[idx])
        except IndexError:
            raise TooHighAmplitudeError, \
                """All the amplitudes are larger than %g. 
                The amplitude position cannot be computed.""" % amplitude

    def relativeAmplitudePosition(self, relAmplitude, timeWindow):
        """Compute the position at which the relative given amplitude is reached.

        The relative amplitude is defined by the amplitude divided by the maximal
        amplitude along the sampled set.
        """
        a = self.fluctuationAmplitudes(timeWindow)
        absAmplitude = relAmplitude*np.max(a)
        return self.amplitudePosition(absAmplitude, timeWindow)

    def maximalAmplitude(self, rMin, rMax, timeWindow):
        """Compute the maximal fluctuation amplitude between rMin and rMax."""
        r = self.positions()
        a = self.fluctuationAmplitudes(timeWindow)
        rMinInt = max(rMin, min(r))
        rMaxInt = min(rMax, max(r))
        f_a = interp1d(r, a, kind='cubic')
        return np.max(f_a(np.linspace(rMinInt, rMaxInt, 100)))

    def amplitudeCylindricalIntegral(self, rMin, rMax, timeWindow):
        """Compute the fluctuation integral between rMin and rMax in cylindrical coordinates.

        The formula is
            \int_{rMin}^{rMax} r*a dr
        """
        r = self.positions()
        a = self.fluctuationAmplitudes(timeWindow)
        rMaxInt = min(rMax, max(r))
        f_ra = interp1d(r, r*a, kind='cubic')
        integral, abserr = quad(f_ra, rMin, rMaxInt)
        return integral

    def integralOscillationRadius(self, rMin, rMax, timeWindow):
        """Integral oscillation radius obtained by integrating between rMin and rMax.

        The formula is
            \sqrt(rMin**2 + 2/maxAmplitude*\int_{rMin}^{rMax} r*a dr)
        """
        a_max = self.maximalAmplitude(rMin, rMax, timeWindow)
        integral = self.amplitudeCylindricalIntegral(rMin, rMax, timeWindow)
        return np.sqrt(rMin**2 + 2*integral/a_max)

    def integralOscillationPenetrationDistance(self, rMin, rMax, timeWindow):
        """Integral oscillation penetration length obtained by integrating between rMin and rMax.

        The formula is
            integralOscillationRadius - rMin

            or

            \sqrt(rMin**2 + 2/maxAmplitude*\int_{rMin}^{rMax} r*a dr) - rMin
        """
        return self.integralOscillationRadius(rMin, rMax, timeWindow) - rMin

    def oscillationIntegral(self, rMin, rMax, timeWindow):
        """Integral of the field oscillation between rMin and rMax.

        The formula is
            2\pi\int_{rMin}^{rMax} r*a dr
        """
        return 2*np.pi*self.amplitudeCylindricalIntegral(rMin, rMax, timeWindow)

    def normalizedOscillationIntegral(self, rMin, rMax, timeWindow):
        """Normalized integral of the field oscillation between rMin and rMax.
        
        The integral is normalized using the maximal oscillation between rMin and rMax as

            2\pi\int_{rMin}^{rMax} r*a dr / (\pi (rMax^2 - rMin^2)*max(a)
        """
        r = self.positions()
        rMaxInt = min(rMax, max(r))
        a_max = self.maximalAmplitude(rMin, rMax, timeWindow)
        return self.oscillationIntegral(rMin, rMax, timeWindow)\
               /(np.pi*(rMaxInt**2 - rMin**2)*a_max)

    def lastMinMax(self, timeWindow=0.0):
        """Compute the minimum and maximum of the sampled set over the time window.

        The extrema are extracted over the interval [lastTime - timeWindow, lastTime].
        If timeWindow equals 0, the last local extrema are computed.
        """
        if timeWindow > 0.0:
            lastTime = self.timeValues()[-1]
            timeKeys = self.timeKeys(lastTime - timeWindow)
            values = np.asarray([self.sampledSetsDict[t] for t in timeKeys])
            return np.min(values, axis=0), np.max(values, axis=0)
        else:
            return self.lastLocalMinMax()

    def lastLocalMinMax(self):
        """Compute the last local minimum and maximum of the sampled set."""
        timeKeys = self.timeKeys()
        values = np.asarray([self.sampledSetsDict[t] for t in timeKeys])
        idxMin = [lastNonStrictLocalExtremum(values[:,i], np.less) 
                                            for i in range(values.shape[1])]
        idxMax = [lastNonStrictLocalExtremum(values[:,i], np.greater) 
                                            for i in range(values.shape[1])]
        minima = [values[idxMin[i], i] for i  in range(values.shape[1])]
        maxima = [values[idxMax[i], i] for i  in range(values.shape[1])]
        return np.asarray(minima), np.asarray(maxima)

    def timeKeys(self, minTime=0.0):
        """
        Return a string list of postprocessing time steps >= minTime.
        """
        times = []
        for key in self.sampledSetsDict:
            try:
                t = float(key)
                if t >= minTime:
                    times.append(key)
            except:
                pass
        return sorted(times) 

    def timeValues(self, minTime=0.0):
        """
        Return a float list of postprocessing time steps >= minTime.
        """
        return [float(t) for t in self.timeKeys(minTime)]


def lastNonStrictLocalExtremum(x, comparator):
    """Return the last non strict local extremum  of the 1D array x."""
    i = firstNonStrictLocalExtremum(x[::-1], comparator)
    return len(x) - i - 1

def firstNonStrictLocalExtremum(x, comparator):
    """Return the first non strict local extremum  of the 1D array x.

    A non-strict extremum is defined as an index i such that there exists 
    positive integers m, n with
        x[i] == x[i-j] for j=1,..,m-1     and  x[i] == x[i+j] for j=1,...,n
        comparator(x[i-m], x[i]) == True  and  comparator(x[i], x[i+n]) == True 

    If no such extremum is found, it returns 0.
    """
    extremumCandidate = False
    iCandidate = 1
    for i in range(1, len(x)):
        if comparator(x[i], x[i-1]):
            extremumCandidate = True
            iCandidate = i
        elif comparator(x[i-1], x[i]):
            if extremumCandidate:
                return iCandidate
            extremumCandidate = False
    return 0


if __name__ == "__main__":
    caseName = '.'
    sampledDirName = 'yProfiles'
    sampledSetFile = 'midstreamProfile_PO2.xy'
    simParams = IOHbO2ParametersAxisymmetric(caseName)
    timeWindow = 2*simParams['RBCLength']/ \
                 (simParams['LDMean']*simParams['RBCVelocity'])
    amplitude = 0.2
    sampledSets = loadSampledSets(caseName, sampledDirName, sampledSetFile,
                                  maxTimeInterval=timeWindow)
    setStats = SampledSetStatistics(sampledSets)
    amplitudes = setStats.fluctuationAmplitudes(timeWindow)
    positions = setStats.sampledSetsDict['positions']
    # print "Position of amplitude = %g mmHg: %g" \
                    # % (amplitude, 
                       # setStats.amplitudePosition(amplitude, timeWindow))
    print "Position of amplitude = %g mmHg: %g" \
                    % (amplitude, 
                       setStats.amplitudePosition(amplitude, timeWindow))
    print 'Relative oscillation radius at 2%: {:g}'. \
        format(setStats.relativeAmplitudePosition(0.02, timeWindow))
    intOscillRad = setStats.integralOscillationRadius(simParams['radiusWall'],
                                                      simParams['radiusTissueLeft'],
                                                      timeWindow)
    print 'Integral oscillation radius: {:g}'.format(intOscillRad)
    # validation code when the oscillation amplitude is replaced by the function 1/r
    # print np.pi*(intOscillRad**2-simParams['radiusWall']**2)/simParams['radiusWall']
    # print 2*np.pi*(max(positions) - simParams['radiusWall'])
    print '{:>12}{:>12}'.format('radius', 'amplitude')
    print '\n'.join(['{:12g}{:12.3g}'.format(r, a) for r, a in zip (positions, amplitudes)])

    minTime = 2.3
    timeKeys = setStats.timeKeys(minTime)
    times = setStats.timeValues(minTime)
    values = np.asarray([setStats.sampledSetsDict[t] for t in timeKeys])
    idxMin = [lastNonStrictLocalExtremum(values[:,i], np.less) 
                                        for i in range(values.shape[1])]
    idxMax = [lastNonStrictLocalExtremum(values[:,i], np.greater) 
                                        for i in range(values.shape[1])]
    minValue, maxValue = setStats.lastMinMax(timeWindow)
    ir = 20
    rPlot = setStats.sampledSetsDict['positions'][ir]
    print "Plotting PO2 fluctuations at r = {:g}".format(rPlot)
    plt.plot(times, values[:,ir])
    plt.plot(times[idxMin[ir]], values[idxMin[ir], ir], 'bo')
    plt.plot(times[idxMax[ir]], values[idxMax[ir], ir], 'ro')
    plt.plot(times[-1], minValue[ir], 'bx')
    plt.plot(times[-1], maxValue[ir], 'rx')

    plt.savefig('plotAmplitudeVsTime_r_{:g}'.format(rPlot))
    plt.show()
    plt.clf()

    print "Plotting PO2 fluctuations as a function of the radius"
    plt.plot(positions, amplitudes)
    plt.plot()
    plt.savefig('plotAmplitudeVsPosition')
    plt.show()

