import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
import argparse
import os
from ast import literal_eval


def plot_mesh(filename, plot_type='pcolor1'):
    parameters = literal_eval(open(filename, 'r').readline().strip())
    xy   = np.loadtxt(filename, delimiter=',', skiprows=1, max_rows=2)
    data = np.loadtxt(filename, delimiter=',', skiprows=3)
    n_x = parameters['n1']
    n_y = parameters['n2']
    x = np.reshape(xy[0, :], (n_x, n_y), order='C')
    y = np.reshape(xy[1, :], (n_x, n_y), order='C')

    status = data[:, 0]

    values = data[:, 1:]

    fig, ax = plt.subplots()


    n_status = len(status[::2])

    def update(i):
        stat = status[2*i]
        val_1 = np.reshape(values[2*i, :], (n_x, n_y), order='C')
        val_2 = np.reshape(values[2*i+1, :], (n_x, n_y), order='C')
        plt.title(stat)
        if plot_type == 'pcolor1':
            plot = plt.pcolor(x, y, val_1, cmap='plasma')
        elif plot_type == 'pcolor2':
            plot = plt.pcolor(x, y, val_2, cmap='plasma')
        elif plot_type == 'vector':
            M = np.hypot(val_1, val_2)
            plot = plt.quiver(x, y, val_1, val_2, M)

        else:
            print("Unkwown plot type!")
            os._exit(1)

        return plot,
    
    anim = FuncAnimation(fig, update, frames=range(n_status), blit=False, repeat=True)
    output = os.path.splitext(filename)[0] + '.gif'
    anim.save(output, writer='imagemagick')


    print("Generated a .gif file at", output)

    """
    H = 20
    R = 1
    r = np.sqrt(x*x + y*y)
    theta = np.arctan2(y, x)
    a = (R-r) * np.pi / H
    phi = np.sin(a) * np.sin(a) * np.cos(2*theta)
    error = abs(val_2 - phi)
    """



parser = argparse.ArgumentParser(description='Projet data processing.')
parser.add_argument('-r', action='store_true', help='run the ./run_project.sh')
parser.add_argument('--plot', help='plot [mesh_u/mesh_v/mesh_w/mesh_p]')
parser.add_argument('--plot_type', help='type of the plot [vector/colorx/colorx] where x is the value to be plotted')

if __name__ == '__main__':
    args = parser.parse_args()
    args = vars(args)
    print(args)
    if args['r']:
        os.system('./run_project.sh')
    if args['plot'] and 'mesh' in args['plot']:
        filename = 'data/' + args['plot'] + '.txt'
        if args['plot_type']:
            plot_mesh(filename, plot_type=args['plot_type'])
        else:
            plot_mesh(filename)

