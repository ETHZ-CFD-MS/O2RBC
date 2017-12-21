#!/usr/bin/env python
#
# Solve the 1D parabolic PDE for PO2 diffusion and consumption in the
# tissue.
#
# Derivation: see CBF, 20.09.2013, "1D Eulerian tissue model"
#
# This script is based on http://jkwiens.com/heat-equation-using-finite-difference/

import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.sparse as sparse
import scipy.sparse.linalg
from   scipy.interpolate import interp1d

from postprocessing import loadProfiles


def load_PO2_pulse():
    # caseName    = '../axisymmetric_RBC/multi_RBC'
    caseName    = '.'
    time        = 0.05
    sampleName  = 'wall'
    u_RBC       = 1.2e-3
    sample_data = loadProfiles.loadProfile(caseName, time, sampleName)

    # transform x-coordinate to time coordinate
    t_values = [(x - sample_data[0][0])/u_RBC for x in sample_data[:,0]]

    # interpolate values
    f_pulse = interp1d(t_values, sample_data[:,1], kind='cubic')

    return f_pulse

def PO2_pulse(f_pulse, t):
    period = max(f_pulse.x) - min(f_pulse.x)
    t_p = t % period
    return f_pulse(t_p)

    # print sample_data

    # the pulse consists in two phases:
    # 1) the pulse itself, that last pulse_fraction*period
    # 2) constant value
    # period = 0.833e-2
    # pulse_fraction = 0.5
    # P_base = 25
    # amplitude = 20
# 
    # t_p = t % period
    # if t_p < pulse_fraction*period:
        # return P_base + amplitude * np.sin(np.pi*t_p/(pulse_fraction*period)) 
    # else:
        # return P_base

def consumption_term(u):
    P_crit = 1
    consumption_rate = 1650
    alpha_tissue = 38.9
    M = consumption_rate / alpha_tissue
    result = np.zeros(len(u))
    for i in xrange(len(u)):
        P = u[i]
        result[i] = -M * P/(P_crit + P)
    return result

def generatePlot(x, u, time):
    # create directory
    dirName = 'plot1DPO2'
    os.system('mkdir -p %s' % dirName)
    i = 0
    plt.plot(1e6*x, u)
    plt.xlabel('r [mum]')
    plt.ylabel('PO2 [mmHg]')
    plt.title('PO2 profile at t = %f' % time)
    plt.ylim([20, 46])
    plt.grid(True)
    plotName = '%s/PO2_t_%.4f.png' % (dirName, time)
    plt.savefig(plotName)
    print 'Saved plot %s.' % plotName
    plt.clf()

def generatePlots(x, data, times):
    # create directory
    dirName = 'plot1DPO2'
    os.system('mkdir -p %s' % dirName)
    i = 0
    for u in data:
        generatePlot(x, u, times[i])
        i = i+1

# parameter
kappa = 2.41e-9
TFinal = 0.2
deltaT = 5e-6
writeInterval = 5e-4

# Domain bounds
R_0 = 3e-6
R_1 = 18e-6

# Load pulse function for boundary condition at inner radius
f_pulse = load_PO2_pulse()

# Initial value
P_init = 25

# Number of internal points
N = 100
 
# Calculate Spatial Step-Size
h = (R_1-R_0)/(N+1.0)
 
# Create number of Time-Steps
NumOfTimeSteps = int(TFinal/deltaT)
writeStepsInterval = int(writeInterval/deltaT)
print 'Total number of time steps: %i' % NumOfTimeSteps
 
# Create grid-points on x axis
x_plot = np.linspace(R_0,R_1,N+2)
x      = x_plot[1:-1]
 
# Initial Conditions
# u = np.transpose(np.mat(10*np.sin(np.pi*x)))
u = P_init*np.ones(N)
 
# Second-Derivative Matrix
data = np.ones((3, N))
data[1] = -2*data[1]
data[1][-1] = -1 # accounts for Neumann boundary condition at exterior radius
diags = [-1,0,1]
D2 = sparse.spdiags(data,diags,N,N)*kappa/(h**2)

# First-Derivative Matrix
data = np.ones((3, N))
data[0] = -1*data[2]
data[1] = 0*data[1]
data[1][-1] = 1 # accounts for Neumann boundary condition at exterior radius
diags = [-1,0,1]
D1 = sparse.spdiags(data,diags,N,N)*kappa/(2*h)
# multiply by 1/x
x_inverse = [1/xx for xx in x]
D_inv = sparse.spdiags(x_inverse, 0, N, N)
D1 = D_inv*D1

D = D2 + D1

# Source term vectors

# consumption term
c = np.zeros(N)
# contribution of second-order derivative
f_a = np.zeros(N)
f_a[0] = PO2_pulse(f_pulse, 0)*kappa/(h**2)

# contribution of first-order derivative
f_b = np.zeros(N)
f_b[0] = -PO2_pulse(f_pulse, 0)*kappa/(2*h*R_0)

# Identity Matrix
I = sparse.identity(N)

# data for each time step
data = [u]
times = [0]
generatePlot(x, u, 0)

for i in range(NumOfTimeSteps):
    t_now = i*deltaT
    t_new = (i+1)*deltaT
    # recompute source term vectors
    f_a[0] =  0.5*(PO2_pulse(f_pulse, t_now)+PO2_pulse(f_pulse, t_new))*kappa/(h**2)
    f_b[0] = -0.5*(PO2_pulse(f_pulse, t_now)+PO2_pulse(f_pulse, t_new))*kappa/(2*h*R_0)
    c      = consumption_term(u)

    # Solve the System: (I - deltaT/2*D2) u_new = (I + deltaT/2*D2)*u_old
    A = (I - deltaT/2*D)
    b = (I + deltaT/2*D)*u + deltaT*(c + f_a + f_b)
    u = sparse.linalg.spsolve( A,  b )

    t_now = t_new

    # Save Data
    if  (i+1) % writeStepsInterval == 0:
        # print 'Saved data at time step %i' % i
        data.append(u)
        times.append(t_now)
        u_plot = np.concatenate(([PO2_pulse(f_pulse, t_now)], u, [u[-1]]))
        generatePlot(x_plot, u_plot, t_now)


