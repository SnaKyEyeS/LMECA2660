import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
import argparse
import os
from ast import literal_eval

def normal(theta, radius):
    U = np.cos(theta) * radius
    V = np.sin(theta) * radius
    return U, V

def tangent(theta, radius):
    U = - np.sin(theta) * radius
    V = np.cos(theta) * radius
    return U, V

def plot_mesh(filename, plot_type='pcolor', value=1, vector=normal):
    parameters = literal_eval(open(filename, 'r').readline().strip())
    xy   = np.loadtxt(filename, delimiter=',', skiprows=1, max_rows=2)
    data = np.loadtxt(filename, delimiter=',', skiprows=3)
    n_x = parameters['n1']
    n_y = parameters['n2']
    x = np.reshape(xy[0, :], (n_x, n_y), order='C')
    y = np.reshape(xy[1, :], (n_x, n_y), order='C')

    theta = np.arctan2(y, x)

    status = data[:, 0]

    values = data[:, 1:]

    fig, _ = plt.subplots()

    plt.xlabel('x')
    plt.ylabel('y')


    n_status = len(status[::2])

    def update(i):
        stat = status[2*i]
        val_1 = np.reshape(values[2*i, :], (n_x, n_y), order='C')
        val_2 = np.reshape(values[2*i+1, :], (n_x, n_y), order='C')

        val = val_1 if value == 1 else val_2

        plt.title(stat)
        if plot_type == 'pcolor':
            plot = plt.pcolor(x, y, val, cmap='plasma')
        elif plot_type == 'vector':
            U, V = vector(theta, val)
            M = np.hypot(U, V)
            plot = plt.quiver(x, y, U, V, M)
        else:
            print("Unkwown plot type!")
            os._exit(1)

        # plt.savefig(f'test{i}.svg')
        return plot,
    
    anim = FuncAnimation(fig, update, frames=range(n_status), blit=False, repeat=True, interval=1000)
    output = os.path.splitext(filename)[0] + '.gif'
    anim.save(output, writer='imagemagick', dpi=300)


    print("Generated a .gif file at", output)



parser = argparse.ArgumentParser(description='Projet data processing.')
parser.add_argument('-r', action='store_true', help='run the ./run_project.sh')
parser.add_argument('--plot', help='plot [mesh_u/mesh_v/mesh_w/mesh_p]')
parser.add_argument('--plot_type', default='pcolor', help='type of the plot [vector/colorx/colorx] where x is the value to be plotted')
parser.add_argument('--value', default='1', help='value to select [1/2]')
parser.add_argument('--vector', default='normal', help='vector field orientation [normal/tangent]')

if __name__ == '__main__':
    args = parser.parse_args()
    args = vars(args)
    print(args)

    filename = 'data/' + args['plot'] + '.txt' if args['plot'] else None
    plot_type = args['plot_type']
    value = int(args['value'])
    vector = normal if args['vector'] == 'normal' else tangent

    if args['r']:
        os.system('./run_project.sh')
    if args['plot']:
        plot_mesh(filename, plot_type=plot_type, value=value, vector=vector)

