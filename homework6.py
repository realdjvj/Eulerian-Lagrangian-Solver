#! /usr/bin/env /Users/vijayn/Code/anaconda2/bin/python

import numpy as np
import csv
import os
import sys
import time
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, FuncFormatter
def plot_centerlines_re(U_100,V_100,U_400,V_400,N,resolution):

    ##### Plotting Function
    X = np.linspace(0.0,1.0,N-1)
    Y = np.linspace(1.0,0.0,N-1)
    globalFontSize = 20

    # set figure and axes properties
    params = {'axes.titlesize'   : globalFontSize,
              'axes.labelsize'   : globalFontSize,
              'figure.autolayout': True,
              'lines.markersize' : 4,
              'lines.linewidth'  : 2.5,
              'xtick.labelsize'  : 15,
              'ytick.labelsize'  : 15,
              'figure.facecolor' : 'white',
              'legend.fontsize'  : globalFontSize,
              'font.weight'      : 'semibold'}
    #          'font.family'      : 'serif',
    #          'font.serif'       : 'Times'}
    plt.rcParams['font.weight'] = 'bold'      # font weight
    plt.rcParams['axes.labelweight'] = 'bold' # axis label weight
    plt.rcParams['axes.linewidth'] = 1.5      # axis linewidth

    # np.seterr('raise')

    plt.figure(figsize=(10, 6))      # figure size
    plt.rc('text', usetex=True)               # use latex things
    plt.rcParams.update(params)               # update parameters
    ax = plt.axes()
    x_minorLocator = matplotlib.ticker.AutoMinorLocator(10) # set tick density for x axis
    y_minorLocator = matplotlib.ticker.AutoMinorLocator(5)  # set tick density for y axis
    ax.xaxis.set_minor_locator(x_minorLocator)
    ax.yaxis.set_minor_locator(y_minorLocator)
    plt.tick_params(axis='both', which='major', direction='in', length=6, width=1,right=True,top=True)
    plt.tick_params(axis='both', which='minor', direction='in', length=3, width=1,right=True,top=True)

    plt.plot(U_100[:,int(N/2)], Y, linewidth=1.5, marker='*', label='Re = 100')
    plt.plot(U_400[:,int(N/2)], Y, linewidth=1.5, marker='^', label='Re = 400')

    ax.set_ylabel(r'$Y/N$',fontsize=globalFontSize)
    ax.set_xlabel(r'$U/U_0$',fontsize=globalFontSize)

    plt.gcf().subplots_adjust(bottom=0.18,left=0.15, right=0.95, top=0.95)

    # plt.xlim(0.0,1.1)
    # plt.ylim(-0.05,0.75)
    # plt.xscale('log')
    lines, labels = ax.get_legend_handles_labels()
    plt.legend(lines, labels, loc='best')
    plt.savefig('./images/x_centerline_%s.png' % (resolution))
    plt.close('all')

    plt.figure(figsize=(10, 6))      # figure size
    plt.rc('text', usetex=True)               # use latex things
    plt.rcParams.update(params)               # update parameters
    ax = plt.axes()
    x_minorLocator = matplotlib.ticker.AutoMinorLocator(10) # set tick density for x axis
    y_minorLocator = matplotlib.ticker.AutoMinorLocator(5)  # set tick density for y axis
    ax.xaxis.set_minor_locator(x_minorLocator)
    ax.yaxis.set_minor_locator(y_minorLocator)
    plt.tick_params(axis='both', which='major', direction='in', length=6, width=1,right=True,top=True)
    plt.tick_params(axis='both', which='minor', direction='in', length=3, width=1,right=True,top=True)

    plt.plot(X, V_100[int(N/2),:], linewidth=1.5, marker='*', label='Re = 100')
    plt.plot(X, V_400[int(N/2),:], linewidth=1.5, marker='*', label='Re = 400')

    ax.set_ylabel(r'$V/U_0$',fontsize=globalFontSize)
    ax.set_xlabel(r'$X/N$',fontsize=globalFontSize)

    plt.gcf().subplots_adjust(bottom=0.18,left=0.15, right=0.95, top=0.95)

    # plt.xlim(0.0,1.1)
    # plt.ylim(-0.05,0.75)
    # plt.xscale('log')
    lines, labels = ax.get_legend_handles_labels()
    plt.legend(lines, labels, loc='best')

    plt.savefig('./images/y_centerline_%s.png' % (resolution))
    plt.close('all')

def plot_centerlines(U,V,N,Re,resolution):

    ##### Plotting Function
    X = np.linspace(0.0,1.0,N-1)
    Y = np.linspace(1.0,0.0,N-1)
    globalFontSize = 20

    # set figure and axes properties
    params = {'axes.titlesize'   : globalFontSize,
              'axes.labelsize'   : globalFontSize,
              'figure.autolayout': True,
              'lines.markersize' : 4,
              'lines.linewidth'  : 2.5,
              'xtick.labelsize'  : 15,
              'ytick.labelsize'  : 15,
              'figure.facecolor' : 'white',
              'legend.fontsize'  : globalFontSize,
              'font.weight'      : 'semibold'}
    #          'font.family'      : 'serif',
    #          'font.serif'       : 'Times'}
    plt.rcParams['font.weight'] = 'bold'      # font weight
    plt.rcParams['axes.labelweight'] = 'bold' # axis label weight
    plt.rcParams['axes.linewidth'] = 1.5      # axis linewidth

    # np.seterr('raise')

    plt.figure(figsize=(10, 6))      # figure size
    plt.rc('text', usetex=True)               # use latex things
    plt.rcParams.update(params)               # update parameters
    ax = plt.axes()
    x_minorLocator = matplotlib.ticker.AutoMinorLocator(10) # set tick density for x axis
    y_minorLocator = matplotlib.ticker.AutoMinorLocator(5)  # set tick density for y axis
    ax.xaxis.set_minor_locator(x_minorLocator)
    ax.yaxis.set_minor_locator(y_minorLocator)
    plt.tick_params(axis='both', which='major', direction='in', length=6, width=1,right=True,top=True)
    plt.tick_params(axis='both', which='minor', direction='in', length=3, width=1,right=True,top=True)

    plt.plot(U[:,int(N/2)], Y, linewidth=1.5, marker='*', label='X-velocity')

    ax.set_ylabel(r'$Y/N$',fontsize=globalFontSize)
    ax.set_xlabel(r'$U/U_0$',fontsize=globalFontSize)

    plt.gcf().subplots_adjust(bottom=0.18,left=0.15, right=0.95, top=0.95)

    # plt.xlim(0.0,1.1)
    # plt.ylim(-0.05,0.75)
    # plt.xscale('log')
    lines, labels = ax.get_legend_handles_labels()
    plt.legend(lines, labels, loc='best')
    plt.savefig('./images/x_centerline_%i_%s.png' % (Re,resolution))
    plt.close('all')

    plt.figure(figsize=(10, 6))      # figure size
    plt.rc('text', usetex=True)               # use latex things
    plt.rcParams.update(params)               # update parameters
    ax = plt.axes()
    x_minorLocator = matplotlib.ticker.AutoMinorLocator(10) # set tick density for x axis
    y_minorLocator = matplotlib.ticker.AutoMinorLocator(5)  # set tick density for y axis
    ax.xaxis.set_minor_locator(x_minorLocator)
    ax.yaxis.set_minor_locator(y_minorLocator)
    plt.tick_params(axis='both', which='major', direction='in', length=6, width=1,right=True,top=True)
    plt.tick_params(axis='both', which='minor', direction='in', length=3, width=1,right=True,top=True)

    plt.plot(X, V[int(N/2),:], linewidth=1.5, marker='*', label='Y-velocity')

    ax.set_ylabel(r'$V/U_0$',fontsize=globalFontSize)
    ax.set_xlabel(r'$X/N$',fontsize=globalFontSize)

    plt.gcf().subplots_adjust(bottom=0.18,left=0.15, right=0.95, top=0.95)

    # plt.xlim(0.0,1.1)
    # plt.ylim(-0.05,0.75)
    # plt.xscale('log')
    lines, labels = ax.get_legend_handles_labels()
    plt.legend(lines, labels, loc='best')

    plt.savefig('./images/y_centerline_%i_%s.png' % (Re,resolution))
    plt.close('all')

def plot_centerlines_3(U_coarse,V_coarse,U_fine,V_fine,U_finer,V_finer,Re,N1,N2,N3):

    ##### Plotting Function
    Y_coarse = np.linspace(1.0,0.0,N1-1)
    Y_fine = np.linspace(1.0,0.0,N2-1)    
    Y_finer = np.linspace(1.0,0.0,N3-1)    
    X_coarse = np.linspace(0.0,1.0,N1-1)
    X_fine   = np.linspace(0.0,1.0,N2-1)
    X_finer  = np.linspace(0.0,1.0,N3-1)
    globalFontSize = 20

    # set figure and axes properties
    params = {'axes.titlesize'   : globalFontSize,
              'axes.labelsize'   : globalFontSize,
              'figure.autolayout': True,
              'lines.markersize' : 4,
              'lines.linewidth'  : 2.5,
              'xtick.labelsize'  : 15,
              'ytick.labelsize'  : 15,
              'figure.facecolor' : 'white',
              'legend.fontsize'  : globalFontSize,
              'font.weight'      : 'semibold'}
    #          'font.family'      : 'serif',
    #          'font.serif'       : 'Times'}
    plt.rcParams['font.weight'] = 'bold'      # font weight
    plt.rcParams['axes.labelweight'] = 'bold' # axis label weight
    plt.rcParams['axes.linewidth'] = 1.5      # axis linewidth

    # np.seterr('raise')

    plt.figure(figsize=(10, 6))      # figure size
    plt.rc('text', usetex=True)               # use latex things
    plt.rcParams.update(params)               # update parameters
    ax = plt.axes()
    x_minorLocator = matplotlib.ticker.AutoMinorLocator(10) # set tick density for x axis
    y_minorLocator = matplotlib.ticker.AutoMinorLocator(5)  # set tick density for y axis
    ax.xaxis.set_minor_locator(x_minorLocator)
    ax.yaxis.set_minor_locator(y_minorLocator)
    plt.tick_params(axis='both', which='major', direction='in', length=6, width=1,right=True,top=True)
    plt.tick_params(axis='both', which='minor', direction='in', length=3, width=1,right=True,top=True)

    plt.plot(U_coarse[:,int(N1/2)], Y_coarse, linewidth=1.5, marker='*', label=r'$\Delta x =$'+'%i'%(N1))
    plt.plot(U_fine[:,int(N2/2)], Y_fine, linewidth=1.5, marker='*',     label=r'$\Delta x =$'+'%i'%(N2))
    plt.plot(U_finer[:,int(N3/2)], Y_finer, linewidth=1.5, marker='*',   label=r'$\Delta x =$'+'%i'%(N3))

    ax.set_ylabel(r'$Y/N$',fontsize=globalFontSize)
    ax.set_xlabel(r'$U/U_0$',fontsize=globalFontSize)

    plt.gcf().subplots_adjust(bottom=0.18,left=0.15, right=0.95, top=0.95)

    # plt.xlim(0.0,1.1)
    # plt.ylim(-0.05,0.75)
    # plt.xscale('log')
    lines, labels = ax.get_legend_handles_labels()
    plt.legend(lines, labels, loc='best')
    plt.savefig('./images/x_centerline_%i_%s.png' % (Re,'3'))
    plt.close('all')

    plt.figure(figsize=(10, 6))      # figure size
    plt.rc('text', usetex=True)               # use latex things
    plt.rcParams.update(params)               # update parameters
    ax = plt.axes()
    x_minorLocator = matplotlib.ticker.AutoMinorLocator(10) # set tick density for x axis
    y_minorLocator = matplotlib.ticker.AutoMinorLocator(5)  # set tick density for y axis
    ax.xaxis.set_minor_locator(x_minorLocator)
    ax.yaxis.set_minor_locator(y_minorLocator)
    plt.tick_params(axis='both', which='major', direction='in', length=6, width=1,right=True,top=True)
    plt.tick_params(axis='both', which='minor', direction='in', length=3, width=1,right=True,top=True)

    plt.plot(X_coarse, V_coarse[int(N1/2),:], linewidth=1.5, marker='*',  label=r'$\Delta y =$'+'%i'%(N1))
    plt.plot(X_fine,   V_fine[int(N2/2),:], linewidth=1.5, marker='*',    label=r'$\Delta y =$'+'%i'%(N2))
    plt.plot(X_finer,  V_finer[int(N3/2),:], linewidth=1.5, marker='*',   label=r'$\Delta y =$'+'%i'%(N3))

    ax.set_ylabel(r'$V/U_0$',fontsize=globalFontSize)
    ax.set_xlabel(r'$X/N$',fontsize=globalFontSize)

    plt.gcf().subplots_adjust(bottom=0.18,left=0.15, right=0.95, top=0.95)

    # plt.xlim(0.0,1.1)
    # plt.ylim(-0.05,0.75)
    # plt.xscale('log')
    lines, labels = ax.get_legend_handles_labels()
    plt.legend(lines, labels, loc='best')

    plt.savefig('./images/y_centerline_%i_%s.png' % (Re,'3'))
    plt.close('all')

def plot_streamlines(U,V,N,Re,resolution):
    ##### Plotting Function
    Y, X = np.mgrid[1.0:0.0:(N-1)*1j, 0.0:1.0:(N-1)*1j]
    globalFontSize = 20
    # set figure and axes properties
    params = {'axes.titlesize'   : globalFontSize,
              'axes.labelsize'   : globalFontSize,
              'figure.autolayout': True,
              'lines.markersize' : 4,
              'lines.linewidth'  : 2.5,
              'xtick.labelsize'  : 15,
              'ytick.labelsize'  : 15,
              'figure.facecolor' : 'white',
              'legend.fontsize'  : globalFontSize,
              'font.weight'      : 'semibold'}
    #          'font.family'      : 'serif',
    #          'font.serif'       : 'Times'}
    plt.rcParams['font.weight'] = 'bold'      # font weight
    plt.rcParams['axes.labelweight'] = 'bold' # axis label weight
    plt.rcParams['axes.linewidth'] = 1.5      # axis linewidth

    # np.seterr('raise')

    plt.figure(figsize=(10, 6))      # figure size
    plt.rc('text', usetex=True)               # use latex things
    plt.rcParams.update(params)               # update parameters
    ax = plt.axes()
    x_minorLocator = matplotlib.ticker.AutoMinorLocator(10) # set tick density for x axis
    y_minorLocator = matplotlib.ticker.AutoMinorLocator(5)  # set tick density for y axis
    ax.xaxis.set_minor_locator(x_minorLocator)
    ax.yaxis.set_minor_locator(y_minorLocator)
    plt.tick_params(axis='both', which='major', direction='in', length=6, width=1,right=True,top=True)
    plt.tick_params(axis='both', which='minor', direction='in', length=3, width=1,right=True,top=True)

    plt.streamplot(X,Y,U,V,density=3.0,linewidth=0.5)

    ax.set_ylabel(r'$Y/N$',fontsize=globalFontSize)
    ax.set_xlabel(r'$X/N$',fontsize=globalFontSize)

    plt.gcf().subplots_adjust(bottom=0.18,left=0.15, right=0.95, top=0.95)

    plt.xlim(0.0,1.0)
    plt.ylim(0.0,1.0)
    # plt.xscale('log')
    # lines, labels = ax.get_legend_handles_labels()
    # plt.legend(lines, labels, loc='best')

    plt.savefig('./images/streamlines_%i_%s.png' % (Re,resolution))
    plt.close('all')

def convergence(U_coarse,U_fine,U_finer):
    error_fine = np.zeros(U_coarse.shape)
    error_coarse = np.zeros(U_coarse.shape)
    error_finer = np.zeros(U_coarse.shape)
    for i in range(0,len(U_coarse[:,0])):
        for j in range(0,len(U_coarse[0,:])):
            error_fine[i,j] = U_finer[4*i,4*j] - U_fine[2*i,2*j]
            error_coarse[i,j] = U_finer[4*i,4*j] - U_coarse[i,j]
            error_finer[i,j] = U_finer[4*i,4*j] - U_finer[4*i,4*j]
    error_fine = np.linalg.norm(error_fine)
    error_coarse = np.linalg.norm(error_coarse)
    error_finer = np.linalg.norm(error_finer) 
    globalFontSize = 20
    # set figure and axes properties
    params = {'axes.titlesize'   : globalFontSize,
              'axes.labelsize'   : globalFontSize,
              'figure.autolayout': True,
              'lines.markersize' : 4,
              'lines.linewidth'  : 2.5,
              'xtick.labelsize'  : 15,
              'ytick.labelsize'  : 15,
              'figure.facecolor' : 'white',
              'legend.fontsize'  : globalFontSize,
              'font.weight'      : 'semibold'}
    #          'font.family'      : 'serif',
    #          'font.serif'       : 'Times'}
    plt.rcParams['font.weight'] = 'bold'      # font weight
    plt.rcParams['axes.labelweight'] = 'bold' # axis label weight
    plt.rcParams['axes.linewidth'] = 1.5      # axis linewidth

    # np.seterr('raise')

    plt.figure(figsize=(10, 6))      # figure size
    plt.rc('text', usetex=True)               # use latex things
    plt.rcParams.update(params)               # update parameters
    ax = plt.axes()
    x_minorLocator = matplotlib.ticker.AutoMinorLocator(10) # set tick density for x axis
    y_minorLocator = matplotlib.ticker.AutoMinorLocator(5)  # set tick density for y axis
    ax.xaxis.set_minor_locator(x_minorLocator)
    ax.yaxis.set_minor_locator(y_minorLocator)
    plt.tick_params(axis='both', which='major', direction='in', length=6, width=1,right=True,top=True)
    plt.tick_params(axis='both', which='minor', direction='in', length=3, width=1,right=True,top=True)

    plt.plot([0.05,0.025,0.0125],[error_coarse,error_fine,error_finer],linewidth=1.5,marker='*')

    ax.set_ylabel(r'Error',fontsize=globalFontSize)
    ax.set_xlabel(r'$\Delta x$',fontsize=globalFontSize)

    plt.gcf().subplots_adjust(bottom=0.18,left=0.15, right=0.95, top=0.95)

    # plt.xlim(1e-1,1e-5)
    # plt.ylim(1e-1,1e-5)
    plt.yscale('log')
    plt.xscale('log')
    # lines, labels = ax.get_legend_handles_labels()
    # plt.legend(lines, labels, loc='best')

    plt.savefig('./images/spatial_%i.png' % Re)
    plt.close('all')

def write_inputs(N,U0,dt,tau):
    with open('input.txt','w') as intxt:
        writer = csv.writer(intxt,delimiter='\t')
        writer.writerow(['N','U0','dt','tau'])
        writer.writerow([N,U0,int(dt),tau])


if __name__ == "__main__":
    if (os.path.exists('./images') != True): os.mkdir('./images') # make images directory if needed
    if (os.path.exists('./data') != True): os.mkdir('./data')

    N1 = 100
    N2 = 2*N1 
    N3 = 2*N2
    U0 = 0.05
    Res = [100, 400, 1000]
    cases = {'coarse':{'N':N1},
             'fine':{'N':N2},
             'finer':{'N':N3}}
    os.system('gfortran -O2 -o lbm.x lbm.f95')
    for case in cases:
        for Re in Res:
            N = cases[case]['N']
            dx = 1.0
            dt = dx**2
            tau = 3.0* U0 * N * dx/Re + 0.5

            write_inputs(N,U0,dt,tau) 
            os.system('./lbm.x')

            U = np.loadtxt('U.csv')
            V = np.loadtxt('V.csv')
            os.system('mv U.csv ./data/U_%i_%s.csv' % (Re,case))
            os.system('mv V.csv ./data/V_%i_%s.csv' % (Re,case))
            plot_streamlines(U,V,N,Re,case)
            plot_centerlines(U,V,N,Re,case)
        U_100 = np.loadtxt('./data/U_100_%s.csv' % case)
        V_100 = np.loadtxt('./data/V_100_%s.csv' % case)
        U_400 = np.loadtxt('./data/U_400_%s.csv' % case)
        V_400 = np.loadtxt('./data/V_400_%s.csv' % case)
        plot_centerlines_re(U_100,V_100,U_400,V_400,N,case)
    for Re in Res:
        U_coarse = np.loadtxt('./data/U_%i_coarse.csv' % Re)
        V_coarse = np.loadtxt('./data/V_%i_coarse.csv' % Re)
        U_fine = np.loadtxt('./data/U_%i_fine.csv' % Re)
        V_fine = np.loadtxt('./data/V_%i_fine.csv' % Re)
        U_finer = np.loadtxt('./data/U_%i_finer.csv' % Re)
        V_finer = np.loadtxt('./data/V_%i_finer.csv' % Re)
        plot_centerlines_3(U_coarse,V_coarse,U_fine,V_fine,U_finer,V_finer,Re,N1,N2,N3)
        convergence(U_coarse,U_fine,U_finer)