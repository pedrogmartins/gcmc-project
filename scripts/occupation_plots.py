import numpy as np
import matplotlib.pyplot as plt
import tqdm
from matplotlib.colors import LogNorm, Normalize
%config InlineBackend.figure_format='retina'


#This scipt extract GCMC trajectories in LAMMPS dump files to
#   create CO2 occupancy 2D histograms in the unit cell.  

%cd trajectories

#Import desired color for histogram
import matplotlib.colors
a = np.asarray(matplotlib.colors.to_rgb('#440154')) - 0.08
matplotlib.colors.to_hex(a)

#Implement plotting function
def plot_2d_occ(file):

    # read the lines in the file
    with open(file) as f:
        lines = f.readlines()

    # Creating a data_frame, columns will be STEP, x, y, and z, fx, fy, fz
    data = np.zeros((len(lines), 4))

    # Parsing through the lines
    index = 0  # Initialize index variable
    for counter, line in enumerate(lines):
        splitline = line.split(' ')

        #Creating a flag that will store the timestep at all times and update as we move down
        if line == 'ITEM: TIMESTEP\n':
            step = int(lines[counter + 1])
            #print(step)

        #Extract x, y coordinates of non-framework CO2 atoms
        if splitline[0].isnumeric() and step < 413850:
            atom_id = int(splitline[0])
            if atom_id > 168 and splitline[2] == 'C':
                data[index, :] = [step, float(splitline[3]), float(splitline[4]),
                                  float(splitline[5])]
                index += 1  # Increment index after adding data

    # Remove unused rows from data
    data = data[:index]

    #Create hetmap for plotting
    steps = data[1:, 0]
    x, y, z = data[1:, 1:4].T
    plt.rcParams['axes.facecolor'] = '#300040'
    heatmap, xedges, yedges = np.histogram2d(x, y, bins = 150)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

    #Plotting Parameters
    fig, ax = plt.subplots(figsize = (20, 20))

    vmin = 0.8  # Minimum value for color mapping
    vmax = 650  # Maximum value for color mapping

    im = ax.imshow(heatmap.T, extent=extent, origin='lower', norm=LogNorm(vmin, vmax))

    # Add colorbar with specified vmin and vmax
    cbar = ax.figure.colorbar(im, ax=ax, shrink=0.4, ticks=[0.8, 1, 10, 100, 650])

    #Add lines around unit cell
    x = np.linspace(0, 100, 1000)
    y = (18.926197347930287)*(x)/(10.927045133565391) + 0.21
    y2  = (18.926197347930287)*(x - 21.854090267130783)/(10.927045133565391) - 0.2

    #im = ax.imshow(heatmap.T, extent=extent, origin='lower', norm=LogNorm())
    #cbar = ax.figure.colorbar(im, ax = ax, shrink=0.4) #, ticks=[0, 0.005, 0.01, 0.015, 0.02])

    for t in cbar.ax.get_yticklabels():
         t.set_fontsize(30)
    plt.xlabel('x (Å)', size = 34)
    plt.ylabel('y (Å)', size = 34)
    plt.plot(x, y, color = "white")
    plt.plot(x, y2, color = "white")
    plt.xticks(fontsize=35)
    plt.yticks(fontsize=35)
    plt.xlim(2.5, 32.5)
    plt.ylim(0, 18)
    fig.tight_layout()
    plt.show()

#Choose range of pressures to plot, this might > 5min if you do all
#   pressures and high occupancy trajectories (loads of CO2 to parse)
ps_selec = [0.001, 0.1, 1, 10]
#[0.001, 0.01, 0.1, 1, 10, 17.7827941, 31.6227766, 56.2341325]

#Run funcion over desired pressure ranges and show plots
for i in ps_selec:
    print(i)
    plot_2d_occ("all_183k_298k_112_p" + str(i) + "_s36825_m5_a20_d0.25_oc0.dat")
