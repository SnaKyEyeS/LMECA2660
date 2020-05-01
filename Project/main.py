import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
import argparse
import os
from ast import literal_eval
import numexpr as ne
from utils.manufactured_solutions import analytical_solutions, parse_solution
import glob

def get_function():
    print("Input the function here.\nIt should take maximum 5 parameters, r (radius), t (theta), R, dt and nu, and return a scalar. Ex.: r*sin(t)\nYou can also use {} around an expression to use a manufactured solution. Ex.: {u_star} + {h_x}\nWrite and press enter to confirm:")
    return parse_solution(input().strip())

def read_file(filename):

    with open(filename, 'r') as file:
        for i, _ in enumerate(file):
            pass
        n_status = (i - 2) >> 1

    file = open(filename, 'r')
    parameters = literal_eval(file.readline().strip())
    xy = np.loadtxt(file, delimiter=',', max_rows=2)

    n_x = parameters['n1']
    n_y = parameters['n2']
    x = np.reshape(xy[0, :], (n_x, n_y), order='C')
    y = np.reshape(xy[1, :], (n_x, n_y), order='C')

    def data_generator(n):
        count = 0
        while count < n:
            data = np.loadtxt(file, delimiter=',', max_rows=2)
            val_1 = data[0, 1:]
            val_2 = data[1, 1:]

            status = data[0, 0]

            yield status, val_1, val_2
            count += 1
        file.close()

    if ('nu' in parameters) and ('dt' in parameters):
        return n_x, n_y, x, y, n_status, data_generator(n_status), parameters['nu'], parameters['dt']
    else:
        return n_x, n_y, x, y, n_status, data_generator(n_status)

def normal(theta, radius):
    U = np.cos(theta) * radius
    V = np.sin(theta) * radius
    return U, V

def tangent(theta, radius):
    U = - np.sin(theta) * radius
    V = np.cos(theta) * radius
    return U, V

def debug_mesh(filename, **kwargs):

    basename = os.path.splitext(os.path.basename(filename))[0]

    n_x, n_y, x, y, _, data_generator, nu, dt = read_file(filename)

    print(f'Loaded nu = {nu} and dt = {dt}')

    analytical = None

    if kwargs['compare']:
        func = get_function()
        r = np.hypot(x, y)
        theta = np.arctan2(y, x)
        analytical = ne.evaluate(func, local_dict={'r':r, 't': theta, 'theta': theta, 'dt': dt, 'nu': nu, 'R': 0.5})

    for status, val_1, val_2 in data_generator:

        val = val_1 if kwargs['value'] == 1 else val_2

        val = np.reshape(val, (n_x, n_y), order='C')

        if kwargs['compare']:
            fig = plt.figure(figsize=plt.figaspect(1/3))
            ax = fig.add_subplot(1, 3, 1, projection='3d')
            ax.plot_surface(x, y, val, cmap='plasma')
            ax.set_title('Results')
            ax = fig.add_subplot(1, 3, 2, projection='3d')
            ax.plot_surface(x, y, analytical, cmap='plasma')
            ax.set_title('Analytical')
            ax = fig.add_subplot(1, 3, 3, projection='3d')
            ax.plot_surface(x, y, np.abs(val-analytical), cmap='plasma')
            ax.set_title('Error')
        else:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

            ax.plot_surface(x, y, val, cmap='plasma')

        plt.suptitle(f'{basename}: t = {status:.5f}s')
        plt.show()

        if not kwargs['debug_all']:
            return

def plot_mesh(filename, **kwargs):

    basename = os.path.splitext(os.path.basename(filename))[0]

    n_x, n_y, x, y, n_status, data_generator = read_file(filename)

    theta = np.arctan2(y, x)

    parse_number = lambda n: str(n).zfill(len(str(n_status)))


    fig, ax = plt.subplots()

    plt.xlabel('x')
    plt.ylabel('y')

    vmin = -5
    vmax =  5

    def set_limits():
        if kwargs['limits']:
            lim = kwargs['limits']
            plt.xlim([-lim, lim])
            plt.ylim([-lim, lim])
        return

    def init():
        set_limits()

    def update(i):

        ax.clear()

        status, val_1, val_2 = next(data_generator)

        val = val_1 if kwargs['value'] == 1 else val_2

        print(f'Generating frame {i+1}/{n_status}')

        global vmin, vmax

        if i == 0:
            vmin = -5
            vmax = 5
            #vmin = 20 * min(val)
            #vmax = 20 * max(val)
            if vmin == vmax:
                vmin = -1
                vmax =  1
        vmin = -5
        vmax = 5

        vmin = vmax = None

        val = np.reshape(val, (n_x, n_y), order='C')

        plt.title(f't = {status:.5f}s')
        if kwargs['plot_type'] == 'pcolor':
            plot = plt.pcolormesh(x, y, val, cmap='Spectral', vmin=vmin, vmax=vmax)
        elif kwargs['plot_type'] == 'vector':
            U, V = kwargs['vector'](theta, val)
            M = np.hypot(U, V)
            plot = plt.quiver(x, y, U, V, M)
        else:
            print("Unkwown plot type!")
            os._exit(1)

        if i == 0:
            fig.colorbar(plot, ax=ax)
        
        plt.autoscale()
        set_limits()

        if kwargs['save_frames']:
            format = kwargs['frame_format']
            plt.savefig(os.path.join(kwargs['output_dir'], f'{basename}_{parse_number(i)}{format}'), dpi=kwargs['dpi'])
        return plot

    anim = FuncAnimation(fig, update, frames=range(n_status), blit=False, init_func=init, repeat=False, interval=kwargs['interval'], cache_frame_data=False)
    
    if kwargs['live']:
        plt.show()
    
    if not kwargs['no_save']:
        if kwargs['movie_format'] == '.mp4':
            output = os.path.join(kwargs['output_dir'], basename + kwargs['movie_format'])
            anim.save(output, writer='ffmpeg', dpi=kwargs['dpi'])
        elif kwargs['movie_format'] == '.gif':
            output = os.path.splitext(filename)[0] + kwargs['movie_format']
            anim.save(output, writer='imagemagick', dpi=kwargs['dpi'])
        else:
            print("Unkwown movie type!")
            os._exit(1)
        print("Generated a {kwargs['movie_format']} file at", output)
    elif kwargs['save_frames']:
        for i in range(n_status):
            update(i)




parser = argparse.ArgumentParser(description='Projet data processing.')
parser.add_argument('-r', action='store_true', help='run the ./run_project.sh')
parser.add_argument('-rebuild', action='store_true', help='tells the ./run_project.sh script to first rebuild and cmake: only needed when new files are added')
parser.add_argument('--plot', help='plot [mesh_u/mesh_v/mesh_w/mesh_p/custom]')
parser.add_argument('--plot_type', default='pcolor', help='type of the plot [vector/pcolor]')
parser.add_argument('--value', default='1', type=int, help='value to select [1/2]')
parser.add_argument('--vector', default='normal', help='vector field orientation [normal/tangent]')
parser.add_argument('-live', action='store_true', help='plot on the fly')
parser.add_argument('-no_save', action='store_true', help='don\'t save the movie of the plot')
parser.add_argument('--interval', default=1000, type=int, help='frame spacing in milliseconds, do not go below 1000 if live option is set!')
parser.add_argument('--dpi', default=300, type=int, help='density per inch, relates to quality')
parser.add_argument('-save_frames', action='store_true', help='save each frame')
parser.add_argument('--movie_format', default='.mp4', help='movie save format [.gif/.mp4]')
parser.add_argument('--frame_format', default='.png', help='frame save format [.png/.svg/...]')
parser.add_argument('--output_dir', default='plots', help='output directory')
parser.add_argument('--input_dir', default='data', help='input directory')
parser.add_argument('--limits', type=float, help='x and y limits, symmetric, positive float or int')
parser.add_argument('-debug', action='store_true', help='enter in debug mode : 3D-plot first status')
parser.add_argument('-debug_all', action='store_true', help='if in debug mode : will 3D-plot all the status')
parser.add_argument('-compare', action='store_true', help='if in debug mode : will ask you to enter the analytical solution to compare with the results')
parser.add_argument('--make_video', help='will make a video from images in .png format located in output directory and named {param}_number.png when {param} is your input')

if __name__ == '__main__':
    args = parser.parse_args()
    args = vars(args)
    print(args)

    filename = args['input_dir'] + '/' + args['plot'] + '.txt' if args['plot'] else None
    args['vector'] = normal if args['vector'] == 'normal' else tangent

    if args['r']:
        if args['rebuild']:
            os.system('./run_project.sh rebuild')
        else:
            os.system('./run_project.sh')
    if args['plot']:
        if args['debug']:
            debug_mesh(filename, **args)
        else:
            plot_mesh(filename, **args)

    if args['make_video']:
        param = args['make_video']
        output_dir = args['output_dir']
        frame_format = args['frame_format']
        filepattern = f'{output_dir}/{param}_*{frame_format}'
        files = glob.glob(filepattern)
        n = files[0][::-1].index('_') - len(frame_format)
        command = f'ffmpeg -r 40 -f image2 -s 1920x1080 -i {output_dir}/{param}_%0{n}d{frame_format} -vcodec libx264 -crf 25  -pix_fmt yuv420p {param}.mp4'
        os.system(command)
        
