import matplotlib.pyplot as plt
import numpy as np

def strouhal_from_diagonostic(filepath, **kwargs):
    data = np.loadtxt(filepath, delimiter=',')

    t_adim = data[:, 0]
    U_inf = 1
    Lc = kwargs['adim_unit_value']
    cd = data[:, 1]
    cl = data[:, 2]
    re_w = data[:, 3]
    y_p = data[:, 4]

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

    ax1.plot(t_adim, cd)
    ax1.set_ylabel('$c_d$')
    ax2.plot(t_adim, cl)
    ax2.set_ylabel('$c_l$')

    adim_symbol = kwargs['adim_unit_symbol']

    xlabel = f'$t U_\infty/' +f'{adim_symbol}$'

    plt.xlabel(xlabel)

    plt.tight_layout()

    plt.savefig('diagnostic_coef.pdf', dpi=kwargs['dpi'])

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

    ax1.plot(t_adim, re_w)
    ax1.set_ylabel('$Re_\omega$')
    ax2.plot(t_adim, cl)
    ax2.set_ylabel('$y^+$')

    adim_symbol = kwargs['adim_unit_symbol']

    xlabel = f'$t U_\infty/' +f'{adim_symbol}$'

    plt.xlabel(xlabel)

    plt.tight_layout()

    plt.savefig('diagnostic_stability.pdf', dpi=kwargs['dpi'])

    index = t_adim > kwargs['time_developed'] # Shedding developed

    t = t_adim[index] * Lc / U_inf

    y = np.fft.fft(cl[index])
    dt = t[1] - t[0]
    n = len(y)
    freq = np.fft.fftfreq(n, d=dt)

    i_max = np.argmax(abs(y))

    f_max = abs(freq[i_max])

    St = f_max * Lc / U_inf

    return f_max, St