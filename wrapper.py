import numpy as np
import subprocess
import csv
import os
import sys
import time
import matplotlib
import matplotlib.pyplot as plt
plt.style.use('tex')
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, FuncFormatter, MaxNLocator

def plot_streamlines(U,V,N,i):
    ##### Plotting Function
    Y, X = np.mgrid[1.0:0.0:(N+1)*1j, 0.0:1.0:(N+1)*1j]

    ax = plt.axes()
    x_minorLocator = matplotlib.ticker.AutoMinorLocator(10) # set tick density for x axis
    y_minorLocator = matplotlib.ticker.AutoMinorLocator(5)  # set tick density for y axis
    ax.xaxis.set_minor_locator(x_minorLocator)
    ax.yaxis.set_minor_locator(y_minorLocator)
    plt.tick_params(axis='both', which='major', direction='in', length=6, width=1,right=True,top=True)
    plt.tick_params(axis='both', which='minor', direction='in', length=3, width=1,right=True,top=True)

    plt.streamplot(X,Y,U,V,density=3.0,linewidth=0.5)

    ax.set_ylabel(r'$Y$')
    ax.set_xlabel(r'$X$')

    plt.gcf().subplots_adjust(bottom=0.18,left=0.15, right=0.95, top=0.95)

    plt.xlim(0.0,1.0)
    plt.ylim(0.0,1.0)
    # plt.xscale('log')
    # lines, labels = ax.get_legend_handles_labels()
    # plt.legend(lines, labels, loc='best')

    plt.savefig('./images/streamlines_{}.pdf'.format(i))
    plt.close('all')

def plot_contours(U,V,N,i):
    ##### Plotting Function
    X = np.linspace(0.0,1.0,N+1)
    Y = np.linspace(1.0,0.0,N+1)

    plt.figure(figsize=(10, 6))      # figure size
    ax = plt.axes()
    x_minorLocator = matplotlib.ticker.AutoMinorLocator(10) # set tick density for x axis
    y_minorLocator = matplotlib.ticker.AutoMinorLocator(5)  # set tick density for y axis
    ax.xaxis.set_minor_locator(x_minorLocator)
    ax.yaxis.set_minor_locator(y_minorLocator)
    plt.tick_params(axis='both', which='major', direction='in', length=6, width=1,right=True,top=True)
    plt.tick_params(axis='both', which='minor', direction='in', length=3, width=1,right=True,top=True)

    # plt.contour(X,Y,U,cmap=)
    fig, ax1 = plt.subplots(nrows=1)
    levels = MaxNLocator(nbins=15).tick_values(-0.1,0.1)
    cf = ax1.contourf(X,Y,U,cmap='coolwarm',vmin=-0.1,vmax=0.1,levels=levels)
    fig.colorbar(cf, ax=ax1)


    ax1.set_ylabel(r'$Y$')
    ax1.set_xlabel(r'$X$')

    # plt.gcf().subplots_adjust(bottom=0.18,left=0.15, right=0.95, top=0.95)

    # plt.xlim(0.0,1.0)
    # plt.ylim(0.0,1.0)
    # plt.xscale('log')
    # lines, labels = ax.get_legend_handles_labels()
    # plt.legend(lines, labels, loc='best')

    plt.savefig('./images/contours_{}.pdf'.format(i))
    plt.close('all')

def plot_scatter(particles,i):
    ##### Plotting Function
    plt.figure(figsize=(10, 6))      # figure size
    ax = plt.axes()
    x_minorLocator = matplotlib.ticker.AutoMinorLocator(10) # set tick density for x axis
    y_minorLocator = matplotlib.ticker.AutoMinorLocator(5)  # set tick density for y axis
    ax.xaxis.set_minor_locator(x_minorLocator)
    ax.yaxis.set_minor_locator(y_minorLocator)
    plt.tick_params(axis='both', which='major', direction='in', length=6, width=1,right=True,top=True)
    plt.tick_params(axis='both', which='minor', direction='in', length=3, width=1,right=True,top=True)

    try:
        plt.scatter(particles[:,1],particles[:,2],s=0.1)
    except:
        print('No more particles')

    ax.set_ylabel(r'$Y$')
    ax.set_xlabel(r'$X$')

    plt.gcf().subplots_adjust(bottom=0.18,left=0.15, right=0.95, top=0.95)

    plt.xlim(0.0,1.0)
    plt.ylim(0.0,1.0)
    # plt.xscale('log')
    # lines, labels = ax.get_legend_handles_labels()
    # plt.legend(lines, labels, loc='best')

    plt.savefig('./images/scatter_{}.pdf'.format(i))
    plt.close('all')

def plot_centerlines(U,V,N,i):
    X = np.linspace(0.0,1.0,N+1)
    Y = np.linspace(1.0,0.0,N+1)
    plt.figure(figsize=(10, 6))      # figure size
    ax = plt.axes()
    x_minorLocator = matplotlib.ticker.AutoMinorLocator(10) # set tick density for x axis
    y_minorLocator = matplotlib.ticker.AutoMinorLocator(5)  # set tick density for y axis
    ax.xaxis.set_minor_locator(x_minorLocator)
    ax.yaxis.set_minor_locator(y_minorLocator)
    plt.tick_params(axis='both', which='major', direction='in', length=6, width=1,right=True,top=True)
    plt.tick_params(axis='both', which='minor', direction='in', length=3, width=1,right=True,top=True)

    plt.plot(U[:,int(N/2)],Y,marker='*',linewidth=1.0)

    ax.set_ylabel(r'$Y$')
    ax.set_xlabel(r'$U$')

    plt.gcf().subplots_adjust(bottom=0.18,left=0.15, right=0.95, top=0.95)

    plt.xlim(-.1,0.1)
    plt.ylim(0.0,1.0)
    # plt.xscale('log')
    # lines, labels = ax.get_legend_handles_labels()
    # plt.legend(lines, labels, loc='best')

    plt.savefig('./images/centerline_x_{}.pdf'.format(i))
    plt.close('all')

    ##### Plotting Function
    plt.figure(figsize=(10, 6))      # figure size
    ax = plt.axes()
    x_minorLocator = matplotlib.ticker.AutoMinorLocator(10) # set tick density for x axis
    y_minorLocator = matplotlib.ticker.AutoMinorLocator(5)  # set tick density for y axis
    ax.xaxis.set_minor_locator(x_minorLocator)
    ax.yaxis.set_minor_locator(y_minorLocator)
    plt.tick_params(axis='both', which='major', direction='in', length=6, width=1,right=True,top=True)
    plt.tick_params(axis='both', which='minor', direction='in', length=3, width=1,right=True,top=True)

    plt.plot(X,U[int(N/2),:],marker='*',linewidth=1.0)
    ax.set_ylabel(r'$U$')
    ax.set_xlabel(r'$X$')

    plt.gcf().subplots_adjust(bottom=0.18,left=0.15, right=0.95, top=0.95)

    plt.xlim(0.0,1.0)
    plt.ylim(-.1,0.1)
    # plt.xscale('log')
    # lines, labels = ax.get_legend_handles_labels()
    # plt.legend(lines, labels, loc='best')

    plt.savefig('./images/centerline_y_{}.pdf'.format(i))
    plt.close('all')

def write_inputs(Re,N,L,U0,dt,jet_height,T_end):
    with open('input.txt','w') as intxt:
        writer = csv.writer(intxt,delimiter='\t')
        writer.writerow(['Re','N','L','U0','dt','jet_height','T_end'])
        writer.writerow([Re,N,L,U0,dt,int(jet_height),T_end])

if __name__ == "__main__":
    if (os.path.exists('./images') != True): os.mkdir('./images') # make images directory if needed
    if (os.path.exists('./data') != True): os.mkdir('./data')
    os.system('gfortran -O2 -ffree-form -c fhdict_module.f90 projection_bak.f95')
    os.system('gfortran -o projection.x projection_bak.o fhdict_module.o')
    L = 1.0
    N = 50
    dx = L/N
    jet_height = int(N/10)
    Re = 1000
    U0 = Re*1.5e-5 / 1.225 /(float(jet_height)/N)
    dt = dx/(10.0*U0)
    # steps = 1000
    T_end = 10.0*(L/U0)
    steps = int(T_end/dt)

    write_inputs(Re,N,L,U0,dt,jet_height,T_end)
    os.system('./projection.x')

    for i in range(0,steps+1,100):
        U = np.loadtxt('./data/U{}.csv'.format(i))
        V = np.loadtxt('./data/V{}.csv'.format(i))
        particles = np.loadtxt('./data/particles{}.csv'.format(i))
        plot_scatter(particles,i)
        plot_streamlines(U,V,N,i)
        plot_contours(U,V,N,i)
        plot_centerlines(U,V,N,i)
        print(i)
    subprocess.run('convert -delay 20 -loop 0 ./images/scatter_{{0..{}..100}}.pdf ./images/scatter.gif'.format(steps), executable='/usr/local/bin/bash', shell=True)
    # subprocess.run('convert -delay 20 -loop 0 ./images/scatter_{{0..{}..100}}.pdf ./images/scatter.gif'.format(steps), executable='/usr/local/bin/bash', shell=True)
    subprocess.run('convert -delay 20 -loop 0 ./images/contours_{{0..{}..100}}.pdf ./images/contours.gif'.format(steps), executable='/usr/local/bin/bash', shell=True)
    subprocess.run('convert -delay 20 -loop 0 ./images/streamlines_{{0..{}..100}}.pdf ./images/streamlines.gif'.format(steps), executable='/usr/local/bin/bash', shell=True)
    subprocess.run('convert -delay 20 -loop 0 ./images/centerline_x_{{0..{}..100}}.pdf ./images/centerline_x.gif'.format(steps), executable='/usr/local/bin/bash', shell=True)
    subprocess.run('convert -delay 20 -loop 0 ./images/centerline_y_{{0..{}..100}}.pdf ./images/centerline_y.gif'.format(steps), executable='/usr/local/bin/bash', shell=True)
