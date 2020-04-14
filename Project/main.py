import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import argparse
import os


def plot_mesh(filename='data/mesh.txt'):
    data = np.loadtxt(filename, delimiter=',')
    n_x = 360
    n_y = 491
    print(data.shape)
    x = np.reshape(data[:, 0], (n_x, n_y))
    y = np.reshape(data[:, 1], (n_x, n_y))
    val_1 = np.reshape(data[:, 2], (n_x, n_y))
    val_2 = np.reshape(data[:, 3], (n_x, n_y))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(x, y, val_2)
    plt.show()



parser = argparse.ArgumentParser(description='Projet data processing.')
parser.add_argument('-r', action='store_true', help='run the ./run_project.sh')
parser.add_argument('--plot', help='plot [mesh/..]')

if __name__ == '__main__':
    args = parser.parse_args()
    args = vars(args)
    print(args)
    if args['r']:
        os.system('./run_project.sh')
    if args['plot'] == 'mesh':
        plot_mesh()

