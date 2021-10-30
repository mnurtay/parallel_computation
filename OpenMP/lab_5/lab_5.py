import numpy as np
import platform
import matplotlib.pyplot as plt
import time


def all_gms(nu, rho, n, x, y, t):
    s = (-1.0) ** (np.ceil(x) + np.ceil(y))

    u = np.pi * s * np.sin(t) * (np.sin(np.pi * x))**2 \
        * np.sin(2.0 * np.pi * y)
    dudt = np.pi * s * np.cos(t) * (np.sin(np.pi * x))**2 \
        * np.sin(2.0 * np.pi * y)
    dudx = np.pi**2 * s * np.sin(t) * np.sin(2.0 * np.pi * x)      \
        * np.sin(2.0 * np.pi * y)
    dudxx = 2.0 * np.pi**3 * s * np.sin(t) * np.cos(2.0 * np.pi * x)      \
        * np.sin(2.0 * np.pi * y)
    dudy = 2.0 * np.pi**2 * s * np.sin(t) * (np.sin(np.pi * x))**2 \
        * np.cos(2.0 * np.pi * y)
    dudyy = - 4.0 * np.pi**3 * s * np.sin(t) * (np.sin(np.pi * x))**2 \
        * np.sin(2.0 * np.pi * y)

    v = - np.pi * s * np.sin(t) * np.sin(2.0 * np.pi * x) \
        * (np.sin(np.pi * y))**2
    dvdt = - np.pi * s * np.cos(t) * np.sin(2.0 * np.pi * x) \
        * (np.sin(np.pi * y))**2
    dvdx = - 2.0 * np.pi**2 * s * np.sin(t) * np.cos(2.0 * np.pi * x) \
        * (np.sin(np.pi * y))**2
    dvdxx = 4.0 * np.pi**3 * s * np.sin(t) * np.sin(2.0 * np.pi * x) \
        * (np.sin(np.pi * y))**2
    dvdy = - np.pi**2 * s * np.sin(t) * np.sin(2.0 * np.pi * x) \
        * np.sin(2.0 * np.pi * y)
    dvdyy = - 2.0 * np.pi**3 * s * np.sin(t) * np.sin(2.0 * np.pi * x) \
        * np.cos(2.0 * np.pi * y)

    p = rho * s * np.sin(t) * np.cos(np.pi * x) \
        * np.sin(np.pi * y)
    dpdt = rho * s * np.cos(t) * np.cos(np.pi * x) \
        * np.sin(np.pi * y)
    dpdx = - np.pi * rho * s * np.sin(t) * np.sin(np.pi * x) \
        * np.sin(np.pi * y)
    dpdxx = - np.pi**2 * rho * s * np.sin(t) * np.cos(np.pi * x) \
        * np.sin(np.pi * y)
    dpdy = np.pi * rho * s * np.sin(t) * np.cos(np.pi * x) \
        * np.cos(np.pi * y)
    dpdyy = - np.pi**2 * rho * s * np.sin(t) * np.cos(np.pi * x) \
        * np.sin(np.pi * y)

    us = dudt + u * dudx + v * dudy + dpdx / rho - nu * (dudxx + dudyy)
    vs = dvdt + u * dvdx + v * dvdy + dpdy / rho - nu * (dvdxx + dvdyy)
    ps = dudx + dvdy

    return \
        u, dudt, dudx, dudxx, dudy, dudyy, us, \
        v, dvdt, dvdx, dvdxx, dvdy, dvdyy, vs, \
        p, dpdt, dpdx, dpdxx, dpdy, dpdyy, ps


def grid_2d(x_num, x_lo, x_hi, y_num, y_lo, y_hi):
    x = np.zeros(x_num * y_num)
    y = np.zeros(x_num * y_num)

    if (x_num == 1):
        xi = (x_lo + x_hi) / 2.0
        k = 0
        for j in range(0, y_num):
            for i in range(0, x_num):
                x[k] = xi
                k = k + 1
    else:
        k = 0
        for j in range(0, y_num):
            for i in range(0, x_num):
                xi = (float(x_num - i - 1) * x_lo
                      + float(i) * x_hi) \
                    / float(x_num - 1)
                x[k] = xi
                k = k + 1

    if (y_num == 1):
        yj = (y_lo + y_hi) / 2.0
        k = 0
        for j in range(0, y_num):
            for i in range(0, x_num):
                y[k] = yj
                k = k + 1
    else:
        k = 0
        for j in range(0, y_num):
            yj = (float(y_num - j - 1) * y_lo
                  + float(j) * y_hi) \
                / float(y_num - 1)
            for i in range(0, x_num):
                y[k] = yj
                k = k + 1

    return x, y


def grid_2d_test():
    x_lo = 10.0
    x_hi = 20.0
    x_num = 5

    y_lo = 5.0
    y_hi = 6.0
    y_num = 3

    [x, y] = grid_2d(x_num, x_lo, x_hi, y_num, y_lo, y_hi)

    k = 0
    for j in range(0, y_num):
        for i in range(0, x_num):

            k = k + 1

    return


def navier_stokes_2d_exact_test():
    r8vec_print_test()

    grid_2d_test()

    uvp_gms_test()
    uvp_gms_test2()
    rhs_gms_test()
    resid_gms_test()
    ns2de_gnuplot_gms_test()
    ns2de_matplotlib_gms_test()

    uvp_lukas_test()
    uvp_lukas_test2()
    rhs_lukas_test()
    resid_lukas_test()
    ns2de_gnuplot_lukas_test()
    ns2de_matplotlib_lukas_test()

    uvp_poiseuille_test()
    uvp_poiseuille_test2()
    rhs_poiseuille_test()
    resid_poiseuille_test()
    ns2de_gnuplot_poiseuille_test()
    ns2de_matplotlib_poiseuille_test()
    parameter_poiseuille_test()

    uvp_spiral_test()
    uvp_spiral_test2()
    rhs_spiral_test()
    resid_spiral_test()
    ns2de_gnuplot_spiral_test()
    ns2de_matplotlib_spiral_test()
    parameter_spiral_test()

    uvp_taylor_test()
    uvp_taylor_test2()
    rhs_taylor_test()
    resid_taylor_test()
    ns2de_gnuplot_taylor_test()
    ns2de_matplotlib_taylor_test()
    parameter_taylor_test()

    print('')
    print('navier_stokes_2d_exact_test:')
    print('  Normal end of execution.')
    return


def ns2de_gnuplot(header, n, x, y, u, v, p, s):
    data_filename = header + '_data.txt'

    data_unit = open(data_filename, 'w')

    for i in range(0, n):
        st = '  %g' % (x[i])
        data_unit.write(st)
        st = '  %g' % (y[i])
        data_unit.write(st)
        st = '  %g' % (u[i])
        data_unit.write(st)
        st = '  %g' % (v[i])
        data_unit.write(st)
        st = '  %g' % (s * u[i])
        data_unit.write(st)
        st = '  %g' % (s * v[i])
        data_unit.write(st)
        st = '  %g' % (p[i])
        data_unit.write(st)
        data_unit.write('\n')

    data_unit.close()

    command_filename = header + '_commands.txt'
    plot_filename = header + '.png'

    command_unit = open(command_filename, 'w')

    command_unit.write('#  %s\n' % (command_filename))
    command_unit.write('#\n')
    command_unit.write('set term png\n')
    command_unit.write('set output "%s"\n' % (plot_filename))
    command_unit.write('#\n')
    command_unit.write('#  Add titles and labels.\n')
    command_unit.write('#\n')
    command_unit.write('set xlabel "<--- X --->"\n')
    command_unit.write('set ylabel "<--- Y --->"\n')
    command_unit.write('set title "Navier-Stokes velocity field"\n')
    command_unit.write('unset key\n')
    command_unit.write('#\n')
    command_unit.write('#  Add grid lines.\n')
    command_unit.write('#\n')
    command_unit.write('set grid\n')
    command_unit.write('set size ratio -1\n')
    command_unit.write('#\n')
    command_unit.write('#  Timestamp the plot.\n')
    command_unit.write('#\n')
    command_unit.write('set timestamp\n')
    command_unit.write(
        'plot "%s" using 1:2:5:6 with vectors \\\n' % (data_filename))
    command_unit.write('  head filled lt 2 linecolor rgb "blue"\n')
    command_unit.write('quit\n')

    data_unit.close()
    return


def ns2de_gnuplot_gms_test():
    nu = 1.0
    rho = 1.0
    t = 1.0

    x_lo = -1.0
    x_hi = 1.0
    x_num = 21

    y_lo = -1.0
    y_hi = 1.0
    y_num = 21

    [x, y] = grid_2d(x_num, x_lo, x_hi, y_num, y_lo, y_hi)

    n = x_num * y_num

    [u, v, p] = uvp_gms(nu, rho, n, x, y, t)

    header = 'gms'
    s = 0.25
    ns2de_gnuplot(header, n, x, y, u, v, p, s)

    return


def ns2de_gnuplot_lukas_test():
    x_lo = 0.0
    x_hi = 1.0
    x_num = 21

    y_lo = 0.0
    y_hi = 1.0
    y_num = 21

    [x, y] = grid_2d(x_num, x_lo, x_hi, y_num, y_lo, y_hi)

    nu = 1.0
    rho = 1.0
    n = x_num * y_num
    t = 0.0

    [u, v, p] = uvp_lukas(nu, rho, n, x, y, t)

    header = 'lukas'
    s = 0.25
    ns2de_gnuplot(header, n, x, y, u, v, p, s)
    return


def ns2de_gnuplot_poiseuille_test():
    x_lo = 0.0
    x_hi = 6.0
    x_num = 61

    y_lo = -1.0
    y_hi = +1.0
    y_num = 21

    [x, y] = grid_2d(x_num, x_lo, x_hi, y_num, y_lo, y_hi)

    nu = 1.0
    rho = 1.0
    n = x_num * y_num
    t = 0.0

    [u, v, p] = uvp_poiseuille(nu, rho, n, x, y, t)

    header = 'poiseuille'
    s = 0.5
    ns2de_gnuplot(header, n, x, y, u, v, p, s)
    return


def ns2de_gnuplot_spiral_test():
    x_lo = 0.0
    x_hi = 1.0
    x_num = 21

    y_lo = 0.0
    y_hi = 1.0
    y_num = 21

    [x, y] = grid_2d(x_num, x_lo, x_hi, y_num, y_lo, y_hi)

    nu = 1.0
    rho = 1.0
    n = x_num * y_num
    t = 0.0

    [u, v, p] = uvp_spiral(nu, rho, n, x, y, t)

    header = 'spiral'
    s = 5.0
    ns2de_gnuplot(header, n, x, y, u, v, p, s)
    return


def ns2de_gnuplot_taylor_test():
    x_lo = 0.5
    x_hi = 2.5
    x_num = 21

    y_lo = 0.5
    y_hi = 2.5
    y_num = 21

    [x, y] = grid_2d(x_num, x_lo, x_hi, y_num, y_lo, y_hi)

    nu = 1.0
    rho = 1.0
    n = x_num * y_num
    t = 0.0

    [u, v, p] = uvp_taylor(nu, rho, n, x, y, t)

    header = 'taylor'
    s = 0.10
    ns2de_gnuplot(header, n, x, y, u, v, p, s)
    return


def ns2de_gnuplot_vortex_test():

    x_lo = 0.5
    x_hi = 1.5
    x_num = 21

    y_lo = 0.5
    y_hi = 1.5
    y_num = 21

    [x, y] = grid_2d(x_num, x_lo, x_hi, y_num, y_lo, y_hi)

    nu = 1.0
    rho = 1.0
    n = x_num * y_num
    t = 0.0

    [u, v, p] = uvp_vortex(nu, rho, n, x, y, t)

    header = 'vortex'
    s = 0.10
    ns2de_gnuplot(header, n, x, y, u, v, p, s)
    return


def ns2de_matplotlib(header, n, x, y, u, v, p, s):
    myplot = plt.figure()
    ax = plt.gca()
    ax.quiver(x, y, u, v)
    ax.set_xlabel('<--X-->')
    ax.set_ylabel('<--Y-->')
    ax.set_title(header)
    ax.axis('equal')
    plt.draw()
    plot_filename = header + '_matplotlib.png'
    myplot.savefig(plot_filename)
    plt.show(block=False)
    plt.close()

    return


def ns2de_matplotlib_gms_test():
    nu = 1.0
    rho = 1.0
    t = 1.0

    x_lo = -1.0
    x_hi = 1.0
    x_num = 21

    y_lo = -1.0
    y_hi = 1.0
    y_num = 21

    x, y = grid_2d(x_num, x_lo, x_hi, y_num, y_lo, y_hi)

    n = x_num * y_num

    u, v, p = uvp_gms(nu, rho, n, x, y, t)

    header = 'gms'
    s = 0.25
    ns2de_matplotlib(header, n, x, y, u, v, p, s)

    return


def ns2de_matplotlib_lukas_test():
    x_lo = 0.0
    x_hi = 1.0
    x_num = 21

    y_lo = 0.0
    y_hi = 1.0
    y_num = 21

    x, y = grid_2d(x_num, x_lo, x_hi, y_num, y_lo, y_hi)

    nu = 1.0
    rho = 1.0
    n = x_num * y_num
    t = 0.0

    u, v, p = uvp_lukas(nu, rho, n, x, y, t)

    header = 'lukas'
    s = 0.25
    ns2de_matplotlib(header, n, x, y, u, v, p, s)
    return


def ns2de_matplotlib_poiseuille_test():
    x_lo = 0.0
    x_hi = 6.0
    x_num = 61

    y_lo = -1.0
    y_hi = +1.0
    y_num = 21

    x, y = grid_2d(x_num, x_lo, x_hi, y_num, y_lo, y_hi)

    nu = 1.0
    rho = 1.0
    n = x_num * y_num
    t = 0.0

    u, v, p = uvp_poiseuille(nu, rho, n, x, y, t)

    header = 'poiseuille'
    s = 5.0
    ns2de_matplotlib(header, n, x, y, u, v, p, s)
    return


def ns2de_matplotlib_spiral_test():
    x_lo = 0.0
    x_hi = 1.0
    x_num = 21

    y_lo = 0.0
    y_hi = 1.0
    y_num = 21

    x, y = grid_2d(x_num, x_lo, x_hi, y_num, y_lo, y_hi)

    nu = 1.0
    rho = 1.0
    n = x_num * y_num
    t = 0.0

    u, v, p = uvp_spiral(nu, rho, n, x, y, t)

    header = 'spiral'
    s = 5.0
    ns2de_matplotlib(header, n, x, y, u, v, p, s)
    return


def ns2de_matplotlib_taylor_test():
    x_lo = 0.5
    x_hi = 2.5
    x_num = 21

    y_lo = 0.5
    y_hi = 2.5
    y_num = 21

    x, y = grid_2d(x_num, x_lo, x_hi, y_num, y_lo, y_hi)

    nu = 1.0
    rho = 1.0
    n = x_num * y_num
    t = 0.0

    u, v, p = uvp_taylor(nu, rho, n, x, y, t)

    header = 'taylor'
    s = 0.10
    ns2de_matplotlib(header, n, x, y, u, v, p, s)
    return


def ns2de_matplotlib_vortex_test():
    x_lo = 0.5
    x_hi = 1.5
    x_num = 21

    y_lo = 0.5
    y_hi = 1.5
    y_num = 21

    x, y = grid_2d(x_num, x_lo, x_hi, y_num, y_lo, y_hi)

    nu = 1.0
    rho = 1.0
    n = x_num * y_num
    t = 0.0

    u, v, p = uvp_vortex(nu, rho, n, x, y, t)

    header = 'vortex'
    s = 0.25
    ns2de_matplotlib(header, n, x, y, u, v, p, s)
    return


def parameter_poiseuille_test():

    print('')
    print('parameter_poiseuille_test')
    print('  Python version: %s' % (platform.python_version()))
    print('  Poiseuille Flow')
    print('  Monitor solution norms for various')
    print('  values of NU, RHO.')

    n = 1000
    x_lo = 0.0
    x_hi = 6.0
    y_lo = -1.0
    y_hi = +1.0
    x = x_lo + (x_hi - x_lo) * np.random.rand(n)
    y = y_lo + (y_hi - y_lo) * np.random.rand(n)
#
#  Vary RHO.
#
    print('')
    print('  RHO affects the pressure scaling.')
    print('')
    print('     RHO         NU           T     ||U||       ||V||       ||P||')
    print('')

    nu = 1.0
    rho = 1.0

    for j in range(0, 3):

        for k in range(0, 6):

            t = k / 5.0

            u, v, p = uvp_poiseuille(nu, rho, n, x, y, t)

            u_norm = np.linalg.norm(u) / n
            v_norm = np.linalg.norm(v) / n
            p_norm = np.linalg.norm(p) / n

            print('  %10.4g  %10.4g  %8.4g  %10.4g  %10.4g  %10.4g'
                  % (rho, nu, t, u_norm, v_norm, p_norm))

        print('')
        rho = rho / 100.0
#
#  Vary NU.
#
    print('')
    print('  NU affects the time scaling.')
    print('')
    print('     RHO         NU           T     ||U||       ||V||       ||P||')
    print('')

    nu = 1.0
    rho = 1.0

    for i in range(0, 4):

        for k in range(0, 6):

            t = k / 5.0

            u, v, p = uvp_poiseuille(nu, rho, n, x, y, t)

            u_norm = np.linalg.norm(u) / n
            v_norm = np.linalg.norm(v) / n
            p_norm = np.linalg.norm(p) / n

            print('  %10.4g  %10.4g  %8.4g  %10.4g  %10.4g  %10.4g'
                  % (rho, nu, t, u_norm, v_norm, p_norm))

        print('')

        nu = nu / 10.0
#
#  Terminate.
#
    print('')
    print('parameter_poiseuille_test:')
    print('  Normal end of execution.')
    return


def parameter_spiral_test():

    print('')
    print('parameter_spiral_test')
    print('  Python version: %s' % (platform.python_version()))
    print('  Spiral Flow')
    print('  Monitor solution norms over time for various')
    print('  values of NU, RHO.')

    n = 1000
    x_lo = 0.0
    x_hi = 1.0
    y_lo = 0.0
    y_hi = 1.0
    x = x_lo + (x_hi - x_lo) * np.random.rand(n)
    y = y_lo + (y_hi - y_lo) * np.random.rand(n)
#
#  Vary RHO.
#
    print('')
    print('  RHO affects the pressure scaling.')
    print('')
    print('     RHO         NU           T     ||U||       ||V||       ||P||')
    print('')

    nu = 1.0
    rho = 1.0

    for j in range(0, 3):

        for k in range(0, 6):

            t = k / 5.0

            u, v, p = uvp_spiral(nu, rho, n, x, y, t)

            u_norm = np.linalg.norm(u) / n
            v_norm = np.linalg.norm(v) / n
            p_norm = np.linalg.norm(p) / n

            print('  %10.4g  %10.4g  %8.4g  %10.4g  %10.4g  %10.4g'
                  % (rho, nu, t, u_norm, v_norm, p_norm))

        print('')
        rho = rho / 100.0
#
#  Vary NU.
#
    print('')
    print('  NU affects the time scaling.')
    print('')
    print('     RHO         NU           T     ||U||       ||V||       ||P||')
    print('')

    nu = 1.0
    rho = 1.0

    for i in range(0, 4):

        for k in range(0, 6):

            t = k / 5.0

            u, v, p = uvp_spiral(nu, rho, n, x, y, t)

            u_norm = np.linalg.norm(u) / n
            v_norm = np.linalg.norm(v) / n
            p_norm = np.linalg.norm(p) / n

            print('  %10.4g  %10.4g  %8.4g  %10.4g  %10.4g  %10.4g'
                  % (rho, nu, t, u_norm, v_norm, p_norm))

        print('')

        nu = nu / 10.0
#
#  Terminate.
#
    print('')
    print('parameter_spiral_test:')
    print('  Normal end of execution.')
    return


def parameter_taylor_test():

    print('')
    print('parameter_taylor_test')
    print('  Python version: %s' % (platform.python_version()))
    print('  Taylor Flow')
    print('  Monitor solution norms over time for various')
    print('  values of NU, RHO.')

    n = 1000
    x_lo = 0.5
    x_hi = 2.5
    y_lo = 0.5
    y_hi = 2.5
    x = x_lo + (x_hi - x_lo) * np.random.rand(n)
    y = y_lo + (y_hi - y_lo) * np.random.rand(n)
#
#  Vary RHO.
#
    print('')
    print('  RHO affects the pressure scaling.')
    print('')
    print('     RHO         NU           T     ||U||       ||V||       ||P||')
    print('')

    nu = 1.0
    rho = 1.0

    for j in range(0, 3):

        for k in range(0, 6):

            t = k / 5.0

            u, v, p = uvp_taylor(nu, rho, n, x, y, t)

            u_norm = np.linalg.norm(u) / n
            v_norm = np.linalg.norm(v) / n
            p_norm = np.linalg.norm(p) / n

            print('  %10.4g  %10.4g  %8.4g  %10.4g  %10.4g  %10.4g'
                  % (rho, nu, t, u_norm, v_norm, p_norm))

        print('')
        rho = rho / 100.0
#
#  Vary NU.
#
    print('')
    print('  NU affects the time scaling.')
    print('')
    print('     RHO         NU           T     ||U||       ||V||       ||P||')
    print('')

    nu = 1.0
    rho = 1.0

    for i in range(0, 4):

        for k in range(0, 6):

            t = k / 5.0

            u, v, p = uvp_taylor(nu, rho, n, x, y, t)

            u_norm = np.linalg.norm(u) / n
            v_norm = np.linalg.norm(v) / n
            p_norm = np.linalg.norm(p) / n

            print('  %10.4g  %10.4g  %8.4g  %10.4g  %10.4g  %10.4g'
                  % (rho, nu, t, u_norm, v_norm, p_norm))

        print('')

        nu = nu / 10.0
#
#  Terminate.
#
    print('')
    print('parameter_taylor_test:')
    print('  Normal end of execution.')
    return


def parameter_vortex_test():

    print('')
    print('parameter_vortex_test')
    print('  Python version: %s' % (platform.python_version()))
    print('  Vortex Flow')
    print('  Monitor solution norms over time for various')
    print('  values of NU, RHO.')

    n = 1000
    x_lo = 0.5
    x_hi = 2.5
    y_lo = 0.5
    y_hi = 2.5
    x = x_lo + (x_hi - x_lo) * np.random.rand(n)
    y = y_lo + (y_hi - y_lo) * np.random.rand(n)
#
#  Vary RHO.
#
    print('')
    print('  RHO affects the pressure scaling.')
    print('')
    print('     RHO         NU           T     ||U||       ||V||       ||P||')
    print('')

    nu = 1.0
    rho = 1.0

    for j in range(0, 3):

        for k in range(0, 6):

            t = k / 5.0

            u, v, p = uvp_vortex(nu, rho, n, x, y, t)

            u_norm = np.linalg.norm(u) / n
            v_norm = np.linalg.norm(v) / n
            p_norm = np.linalg.norm(p) / n

            print('  %10.4g  %10.4g  %8.4g  %10.4g  %10.4g  %10.4g'
                  % (rho, nu, t, u_norm, v_norm, p_norm))

        print('')
        rho = rho / 100.0
#
#  Vary NU.
#
    print('')
    print('  NU affects the time scaling.')
    print('')
    print('     RHO         NU           T     ||U||       ||V||       ||P||')
    print('')

    nu = 1.0
    rho = 1.0

    for i in range(0, 4):

        for k in range(0, 6):

            t = k / 5.0

            u, v, p = uvp_vortex(nu, rho, n, x, y, t)

            u_norm = np.linalg.norm(u) / n
            v_norm = np.linalg.norm(v) / n
            p_norm = np.linalg.norm(p) / n

            print('  %10.4g  %10.4g  %8.4g  %10.4g  %10.4g  %10.4g'
                  % (rho, nu, t, u_norm, v_norm, p_norm))

        print('')

        nu = nu / 10.0
#
#  Terminate.
#
    print('')
    print('parameter_vortex_test:')
    print('  Normal end of execution.')
    return


def r8vec_print(n, a, title):

    print('')
    print(title)
    print('')
    for i in range(0, n):
        print('%6d:  %12g' % (i, a[i]))


def r8vec_print_test():

    print('')
    print('r8vec_print_test')
    print('  Python version: %s' % (platform.python_version()))
    print('  r8vec_print prints an R8VEC.')

    n = 4
    v = np.array([123.456, 0.000005, -1.0E+06, 3.14159265], dtype=np.float64)
    r8vec_print(n, v, '  Here is an R8VEC:')
#
#  Terminate.
#
    print('')
    print('r8vec_print_test:')
    print('  Normal end of execution.')
    return


def resid_gms(nu, rho, n, x, y, t):

    u, dudt, dudx, dudxx, dudy, dudyy, us, \
        v, dvdt, dvdx, dvdxx, dvdy, dvdyy, vs, \
        p, dpdt, dpdx, dpdxx, dpdy, dpdyy, ps = all_gms(nu, rho, n, x, y, t)

    ur = dudt + u * dudx + v * dudy + dpdx / rho - nu * (dudxx + dudyy) - us
    vr = dvdt + u * dvdx + v * dvdy + dpdy / rho - nu * (dvdxx + dvdyy) - vs
    pr = dudx + dvdy - ps

    return ur, vr, pr


def resid_gms_test():

    nu = 1.0
    rho = 1.0
    t = 1.0

    print('')
    print('resid_gms_test')
    print('  Python version: %s' % (platform.python_version()))
    print('  GMS Flow')
    print('  Sample the Navier-Stokes residuals')
    print('  over the [-1,+1]x[-1,+1] square')
    print('  at time T = 1.0')
    print('  Kinematic viscosity NU = %g' % (nu))
    print('  Fluid density RHO = %g' % (rho))

    n = 1000
    x_lo = -1.0
    x_hi = 1.0
    y_lo = -1.0
    y_hi = 1.0
    x = x_lo + (x_hi - x_lo) * np.random.rand(n)
    y = y_lo + (y_hi - y_lo) * np.random.rand(n)

    ur, vr, pr = resid_gms(nu, rho, n, x, y, t)

    print('')
    print('           Minimum       Maximum')
    print('')
    print('  Ur:  %14.6g  %14.6g' % (np.min(np.abs(ur)), np.max(np.abs(ur))))
    print('  Vr:  %14.6g  %14.6g' % (np.min(np.abs(vr)), np.max(np.abs(vr))))
    print('  Pr:  %14.6g  %14.6g' % (np.min(np.abs(pr)), np.max(np.abs(pr))))
#
#  Terminate.
#
    print('')
    print('resid_gms_test:')
    print('  Normal end of execution.')
    return


def resid_lukas(nu, rho, n, x, y, t):

    ur = np.zeros(n)
    vr = np.zeros(n)
    pr = np.zeros(n)
#
#  Get the right hand side functions.
#
    f, g, h = rhs_lukas(nu, rho, n, x, y, t)
#
#  Form the functions and derivatives for the left hand side.
#
    u = - np.cos(np.pi * x) / np.pi
    dudt = np.zeros(n)
    dudx = np.sin(np.pi * x)
    dudxx = np.pi * np.cos(np.pi * x)
    dudy = np.zeros(n)
    dudyy = np.zeros(n)

    v = - y * np.sin(np.pi * x)
    dvdt = np.zeros(n)
    dvdx = - np.pi * y * np.cos(np.pi * x)
    dvdxx = + np.pi * np.pi * y * np.sin(np.pi * x)
    dvdy = - np.sin(np.pi * x)
    dvdyy = np.zeros(n)

    p = np.zeros(n)
    dpdx = np.zeros(n)
    dpdy = np.zeros(n)
#
#  Evaluate the residuals.
#
    ur = dudt - nu * (dudxx + dudyy) + u * dudx + v * dudy + dpdx / rho - f
    vr = dvdt - nu * (dvdxx + dvdyy) + u * dvdx + v * dvdy + dpdy / rho - g
    pr = dudx + dvdy - h

    return ur, vr, pr


def resid_lukas_test():

    nu = 1.0
    rho = 1.0

    print('')
    print('resid_lukas_test')
    print('  Python version: %s' % (platform.python_version()))
    print('  Lukas Bystricky Flow')
    print('  Sample the Navier-Stokes residuals')
    print('  at the initial time T = 0, over the unit square.')
    print('  Kinematic viscosity NU = %g' % (nu))
    print('  Fluid density RHO = %g' % (rho))

    n = 1000
    x_lo = 0.0
    x_hi = 1.0
    y_lo = 0.0
    y_hi = 1.0
    x = x_lo + (x_hi - x_lo) * np.random.rand(n)
    y = y_lo + (y_hi - y_lo) * np.random.rand(n)
    t = 0.0

    ur, vr, pr = resid_lukas(nu, rho, n, x, y, t)

    print('')
    print('           Minimum       Maximum')
    print('')
    print('  Ur:  %14.6g  %14.6g' % (np.min(np.abs(ur)), np.max(np.abs(ur))))
    print('  Vr:  %14.6g  %14.6g' % (np.min(np.abs(vr)), np.max(np.abs(vr))))
    print('  Pr:  %14.6g  %14.6g' % (np.min(np.abs(pr)), np.max(np.abs(pr))))
#
#  Terminate.
#
    print('')
    print('resid_lukas_test:')
    print('  Normal end of execution.')
    return


def resid_poiseuille(nu, rho, n, x, y, t):

    ur = np.zeros(n)
    vr = np.zeros(n)
    pr = np.zeros(n)
#
#  Get the right hand side functions.
#
    f, g, h = rhs_poiseuille(nu, rho, n, x, y, t)
#
#  Form the functions and derivatives for the left hand side.
#
    u = 1.0 - y ** 2
    dudt = 0.0
    dudx = 0.0
    dudxx = 0.0
    dudy = - 2.0 * y
    dudyy = - 2.0

    v = 0.0
    dvdt = 0.0
    dvdx = 0.0
    dvdxx = 0.0
    dvdy = 0.0
    dvdyy = 0.0

    p = - 2.0 * nu * rho * x
    dpdx = - 2.0 * nu * rho
    dpdy = 0.0
#
#  Evaluate the residuals.
#
    ur = dudt - nu * (dudxx + dudyy) \
        + u * dudx + v * dudy + dpdx / rho - f

    vr = dvdt - nu * (dvdxx + dvdyy) \
        + u * dvdx + v * dvdy + dpdy / rho - g

    pr = dudx + dvdy - h

    return ur, vr, pr


def resid_poiseuille_test():

    nu = 1.0
    rho = 1.0

    print('')
    print('resid_poiseuille_test')
    print('  Python version: %s' % (platform.python_version()))
    print('  Poiseuille Flow:')
    print('  Sample the Navier-Stokes residuals')
    print('  at the initial time T = 0, over a channel region.')
    print('  Kinematic viscosity NU = %g' % (nu))
    print('  Fluid density RHO = %g' % (rho))

    n = 1000
    x_lo = 0.0
    x_hi = 6.0
    y_lo = -1.0
    y_hi = +1.0
    x = x_lo + (x_hi - x_lo) * np.random.rand(n)
    y = y_lo + (y_hi - y_lo) * np.random.rand(n)
    t = 0.0

    ur, vr, pr = resid_poiseuille(nu, rho, n, x, y, t)

    print('')
    print('           Minimum       Maximum')
    print('')
    print('  Ur:  %14.6g  %14.6g' % (np.min(np.abs(ur)), np.max(np.abs(ur))))
    print('  Vr:  %14.6g  %14.6g' % (np.min(np.abs(vr)), np.max(np.abs(vr))))
    print('  Pr:  %14.6g  %14.6g' % (np.min(np.abs(pr)), np.max(np.abs(pr))))
#
#  Terminate.
#
    print('')
    print('resid_poiseuille_test:')
    print('  Normal end of execution.')
    return


def resid_spiral(nu, rho, n, x, y, t):

    ur = np.zeros(n)
    vr = np.zeros(n)
    pr = np.zeros(n)
#
#  Get the right hand side functions.
#
    f, g, h = rhs_spiral(nu, rho, n, x, y, t)
#
#  Form the functions and derivatives for the left hand side.
#
    u = (1.0 + nu * t) * 2.0 \
        * (x ** 4 - 2.0 * x ** 3 + x ** 2) \
        * (2.0 * y ** 3 - 3.0 * y ** 2 + y)

    dudt = nu * 2.0 \
        * (x ** 4 - 2.0 * x ** 3 + x ** 2) \
        * (2.0 * y ** 3 - 3.0 * y ** 2 + y)

    dudx = (1.0 + nu * t) * 2.0 \
        * (4.0 * x ** 3 - 6.0 * x ** 2 + 2.0 * x) \
        * (2.0 * y ** 3 - 3.0 * y ** 2 + y)

    dudxx = (1.0 + nu * t) * 2.0 \
        * (12.0 * x ** 2 - 12.0 * x + 2.0) \
        * (2.0 * y ** 3 - 3.0 * y ** 2 + y)

    dudy = (1.0 + nu * t) * 2.0 \
        * (x ** 4 - 2.0 * x ** 3 + x ** 2) \
        * (6.0 * y ** 2 - 6.0 * y + 1.0)

    dudyy = (1.0 + nu * t) * 2.0 \
        * (x ** 4 - 2.0 * x ** 3 + x ** 2) \
        * (12.0 * y - 6.0)

    v = - (1.0 + nu * t) * 2.0 \
        * (2.0 * x ** 3 - 3.0 * x ** 2 + x) \
        * (y ** 4 - 2.0 * y ** 3 + y ** 2)

    dvdt = - nu * 2.0 \
        * (2.0 * x ** 3 - 3.0 * x ** 2 + x) \
        * (y ** 4 - 2.0 * y ** 3 + y ** 2)

    dvdx = - (1.0 + nu * t) * 2.0 \
        * (6.0 * x ** 2 - 6.0 * x + 1.0) \
        * (y ** 4 - 2.0 * y ** 3 + y ** 2)

    dvdxx = - (1.0 + nu * t) * 2.0 \
        * (12.0 * x - 6.0) \
        * (y ** 4 - 2.0 * y ** 3 + y ** 2)

    dvdy = - (1.0 + nu * t) * 2.0 \
        * (2.0 * x ** 3 - 3.0 * x ** 2 + x) \
        * (4.0 * y ** 3 - 6.0 * y ** 2 + 2.0 * y)

    dvdyy = - (1.0 + nu * t) * 2.0 \
        * (2.0 * x ** 3 - 3.0 * x ** 2 + x) \
        * (12.0 * y ** 2 - 12.0 * y + 2.0)

    p = rho * y
    dpdx = 0.0
    dpdy = rho
#
#  Evaluate the residuals.
#
    ur = dudt - nu * (dudxx + dudyy) \
        + u * dudx + v * dudy + dpdx / rho - f

    vr = dvdt - nu * (dvdxx + dvdyy) \
        + u * dvdx + v * dvdy + dpdy / rho - g

    pr = dudx + dvdy - h

    return ur, vr, pr


def resid_spiral_test():

    nu = 1.0
    rho = 1.0

    print('')
    print('resid_spiral_test')
    print('  Python version: %s' % (platform.python_version()))
    print('  Spiral Flow:')
    print('  Sample the Navier-Stokes residuals')
    print('  at the initial time T = 0, over the unit square.')
    print('  Kinematic viscosity NU = %g' % (nu))
    print('  Fluid density RHO = %g' % (rho))

    n = 1000
    x_lo = 0.0
    x_hi = 1.0
    y_lo = 0.0
    y_hi = 1.0
    x = x_lo + (x_hi - x_lo) * np.random.rand(n)
    y = y_lo + (y_hi - y_lo) * np.random.rand(n)
    t = 0.0

    ur, vr, pr = resid_spiral(nu, rho, n, x, y, t)

    print('')
    print('           Minimum       Maximum')
    print('')
    print('  Ur:  %14.6g  %14.6g' % (np.min(np.abs(ur)), np.max(np.abs(ur))))
    print('  Vr:  %14.6g  %14.6g' % (np.min(np.abs(vr)), np.max(np.abs(vr))))
    print('  Pr:  %14.6g  %14.6g' % (np.min(np.abs(pr)), np.max(np.abs(pr))))
#
#  Terminate.
#
    print('')
    print('resid_spiral_test:')
    print('  Normal end of execution.')
    return


def resid_taylor(nu, rho, n, x, y, t):

    #
    #  Get the right hand sides.
    #
    f, g, h = rhs_taylor(nu, rho, n, x, y, t)
#
#  Make space.
#
    c2x = np.array(n)
    c2y = np.array(n)
    cx = np.array(n)
    cy = np.array(n)
    e2t = np.array(n)
    e4t = np.array(n)
    p = np.array(n)
    px = np.array(n)
    py = np.array(n)
    s2x = np.array(n)
    s2y = np.array(n)
    sx = np.array(n)
    sy = np.array(n)
    u = np.array(n)
    ut = np.array(n)
    ux = np.array(n)
    uxx = np.array(n)
    uy = np.array(n)
    uyy = np.array(n)
    v = np.array(n)
    vt = np.array(n)
    vx = np.array(n)
    vxx = np.array(n)
    vy = np.array(n)
    vyy = np.array(n)
#
#  Make some temporaries.
#
    cx = np.cos(np.pi * x)
    cy = np.cos(np.pi * y)

    sx = np.sin(np.pi * x)
    sy = np.sin(np.pi * y)

    e2t = np.exp(- 2.0 * np.pi * np.pi * nu * t)

    c2x = np.cos(2.0 * np.pi * x)
    c2y = np.cos(2.0 * np.pi * y)

    s2x = np.sin(2.0 * np.pi * x)
    s2y = np.sin(2.0 * np.pi * y)

    e4t = np.exp(- 4.0 * np.pi * np.pi * nu * t)
#
#  Form the functions and derivatives.
#
    u = -                            cx * sy * e2t
    dudx = np.pi * sx * sy * e2t
    dudxx = np.pi * np.pi * cx * sy * e2t
    dudy = -                    np.pi * cx * cy * e2t
    dudyy = np.pi * np.pi * cx * sy * e2t
    dudt = + 2.0 * nu * np.pi * np.pi * cx * sy * e2t

    v = sx * cy * e2t
    dvdx = np.pi * cx * cy * e2t
    dvdxx = -            np.pi * np.pi * sx * cy * e2t
    dvdy = -                    np.pi * sx * sy * e2t
    dvdyy = -            np.pi * np.pi * sx * cy * e2t
    dvdt = - 2.0 * nu * np.pi * np.pi * sx * cy * e2t

    p = - 0.25 * rho * (c2x + c2y) * e4t
    dpdx = + 0.5 * rho * np.pi * s2x * e4t
    dpdy = + 0.5 * rho * np.pi * s2y * e4t
#
#  Evaluate the residuals.
#
    ur = dudt + u * dudx + v * dudy + (1.0 / rho) * dpdx \
        - nu * (dudxx + dudyy) - f

    vr = dvdt + u * dvdx + v * dvdy + (1.0 / rho) * dpdy \
        - nu * (dvdxx + dvdyy) - g

    pr = dudx + dvdy - h

    return ur, vr, pr


def resid_taylor_test():

    nu = 1.0
    rho = 1.0

    print('')
    print('resid_taylor_test')
    print('  Python version: %s' % (platform.python_version()))
    print('  Sample the Taylor residuals')
    print('  at the initial time T = 0, using a region that is')
    print('  the square centered at (1.5,1.5) with "radius" 1.0,')
    print('  Kinematic viscosity NU = %g' % (nu))
    print('  Fluid density RHO = %g' % (rho))

    n = 1000
    x_lo = 0.5
    x_hi = +2.5
    y_lo = 0.5
    y_hi = 2.5
    x = x_lo + (x_hi - x_lo) * np.random.rand(n)
    y = y_lo + (y_hi - y_lo) * np.random.rand(n)
    t = 0.0

    ur, vr, pr = resid_taylor(nu, rho, n, x, y, t)

    print('')
    print('           Minimum       Maximum')
    print('')
    print('  Ur:  %14.6g  %14.6g' % (np.min(np.abs(ur)), np.max(np.abs(ur))))
    print('  Vr:  %14.6g  %14.6g' % (np.min(np.abs(vr)), np.max(np.abs(vr))))
    print('  Pr:  %14.6g  %14.6g' % (np.min(np.abs(pr)), np.max(np.abs(pr))))
#
#  Terminate.
#
    print('')
    print('resid_taylor_test:')
    print('  Normal end of execution.')
    return


def resid_vortex(nu, rho, n, x, y, t):

    #
    #  Get the right hand sides.
    #
    f, g, h = rhs_vortex(nu, rho, n, x, y, t)
#
#  Make space.
#
    c2x = np.array(n)
    c2y = np.array(n)
    cx = np.array(n)
    cy = np.array(n)
    e2t = np.array(n)
    e4t = np.array(n)
    p = np.array(n)
    px = np.array(n)
    py = np.array(n)
    s2x = np.array(n)
    s2y = np.array(n)
    sx = np.array(n)
    sy = np.array(n)
    u = np.array(n)
    ut = np.array(n)
    ux = np.array(n)
    uxx = np.array(n)
    uy = np.array(n)
    uyy = np.array(n)
    v = np.array(n)
    vt = np.array(n)
    vx = np.array(n)
    vxx = np.array(n)
    vy = np.array(n)
    vyy = np.array(n)
#
#  Make some temporaries.
#
    cx = np.cos(np.pi * x)
    cy = np.cos(np.pi * y)

    sx = np.sin(np.pi * x)
    sy = np.sin(np.pi * y)

    c2x = np.cos(2.0 * np.pi * x)
    c2y = np.cos(2.0 * np.pi * y)

    s2x = np.sin(2.0 * np.pi * x)
    s2y = np.sin(2.0 * np.pi * y)
#
#  Form the functions and derivatives.
#
    u = -                            cx * sy
    dudx = np.pi * sx * sy
    dudxx = np.pi * np.pi * cx * sy
    dudy = -                    np.pi * cx * cy
    dudyy = np.pi * np.pi * cx * sy
    dudt = np.zeros(n)

    v = sx * cy
    dvdx = np.pi * cx * cy
    dvdxx = -            np.pi * np.pi * sx * cy
    dvdy = -                    np.pi * sx * sy
    dvdyy = -            np.pi * np.pi * sx * cy
    dvdt = np.zeros(n)

    p = - 0.25 * rho * (c2x + c2y)
    dpdx = + 0.5 * rho * np.pi * s2x
    dpdy = + 0.5 * rho * np.pi * s2y
#
#  Evaluate the residuals.
#
    ur = dudt + u * dudx + v * dudy + (1.0 / rho) * dpdx \
        - nu * (dudxx + dudyy) - f

    vr = dvdt + u * dvdx + v * dvdy + (1.0 / rho) * dpdy \
        - nu * (dvdxx + dvdyy) - g

    pr = dudx + dvdy - h

    return ur, vr, pr


def resid_vortex_test():

    nu = 1.0
    rho = 1.0

    print('')
    print('resid_vortex_test')
    print('  Python version: %s' % (platform.python_version()))
    print('  Sample the Vortex residuals')
    print('  at the initial time T = 0, using a region that is')
    print('  the square centered at (1.5,1.5) with "radius" 1.0,')
    print('  Kinematic viscosity NU = %g' % (nu))
    print('  Fluid density RHO = %g' % (rho))

    n = 1000
    x_lo = 0.5
    x_hi = +2.5
    y_lo = 0.5
    y_hi = 2.5
    x = x_lo + (x_hi - x_lo) * np.random.rand(n)
    y = y_lo + (y_hi - y_lo) * np.random.rand(n)
    t = 0.0

    ur, vr, pr = resid_vortex(nu, rho, n, x, y, t)

    print('')
    print('           Minimum       Maximum')
    print('')
    print('  Ur:  %14.6g  %14.6g' % (np.min(np.abs(ur)), np.max(np.abs(ur))))
    print('  Vr:  %14.6g  %14.6g' % (np.min(np.abs(vr)), np.max(np.abs(vr))))
    print('  Pr:  %14.6g  %14.6g' % (np.min(np.abs(pr)), np.max(np.abs(pr))))
#
#  Terminate.
#
    print('')
    print('resid_vortex_test:')
    print('  Normal end of execution.')
    return


def rhs_gms(nu, rho, n, x, y, t):

    u, dudt, dudx, dudxx, dudy, dudyy, us, \
        v, dvdt, dvdx, dvdxx, dvdy, dvdyy, vs, \
        p, dpdt, dpdx, dpdxx, dpdy, dpdyy, ps = all_gms(nu, rho, n, x, y, t)

    return us, vs, ps


def rhs_gms_test():

    nu = 1.0
    rho = 1.0
    t = 1.0

    print('')
    print('rhs_gms_test')
    print('  Python version: %s' % (platform.python_version()))
    print('  GMS Flow')
    print('  Sample the Navier-Stokes right hand sides')
    print('  over the [-1,+1]x[-1,+1] square')
    print('  at time T = 1.0')
    print('  Kinematic viscosity NU = %g' % (nu))
    print('  Fluid density RHO = %g' % (rho))

    n = 1000
    x_lo = -1.0
    x_hi = 1.0
    y_lo = -1.0
    y_hi = 1.0
    x = x_lo + (x_hi - x_lo) * np.random.rand(n)
    y = y_lo + (y_hi - y_lo) * np.random.rand(n)
    t = 0.0

    f, g, h = rhs_gms(nu, rho, n, x, y, t)

    print('')
    print('           Minimum       Maximum')
    print('')
    print('  Ur:  %14.6g  %14.6g' % (np.min(f), np.max(f)))
    print('  Vr:  %14.6g  %14.6g' % (np.min(g), np.max(g)))
    print('  Pr:  %14.6g  %14.6g' % (np.min(h), np.max(h)))
#
#  Terminate.
#
    print('')
    print('rhs_gms_test:')
    print('  Normal end of execution.')
    return


def rhs_lukas(nu, rho, n, x, y, t):

    f = np.zeros(n)
    g = np.zeros(n)
    h = np.zeros(n)

    u = - np.cos(np.pi * x) / np.pi
    dudt = np.zeros(n)
    dudx = np.sin(np.pi * x)
    dudxx = np.pi * np.cos(np.pi * x)
    dudy = np.zeros(n)
    dudyy = np.zeros(n)

    v = - y * np.sin(np.pi * x)
    dvdt = np.zeros(n)
    dvdx = - np.pi * y * np.cos(np.pi * x)
    dvdxx = + np.pi * np.pi * y * np.sin(np.pi * x)
    dvdy = - np.sin(np.pi * x)
    dvdyy = np.zeros(n)

    p = np.zeros(n)
    dpdx = np.zeros(n)
    dpdy = np.zeros(n)

    f = dudt - nu * (dudxx + dudyy) + u * dudx + v * dudy + dpdx / rho
    g = dvdt - nu * (dvdxx + dvdyy) + u * dvdx + v * dvdy + dpdy / rho
    h = dudx + dvdy

    return f, g, h


def rhs_lukas_test():

    nu = 1.0
    rho = 1.0

    print('')
    print('rhs_lukas_test')
    print('  Python version: %s' % (platform.python_version()))
    print('  Lukas Bystricky Flow')
    print('  Sample the Navier-Stokes right hand sides')
    print('  at the initial time T = 0, over the unit square.')
    print('  Kinematic viscosity NU = %g' % (nu))
    print('  Fluid density RHO = %g' % (rho))

    n = 1000
    x_lo = 0.0
    x_hi = 1.0
    y_lo = 0.0
    y_hi = 1.0
    x = x_lo + (x_hi - x_lo) * np.random.rand(n)
    y = y_lo + (y_hi - y_lo) * np.random.rand(n)
    t = 0.0

    f, g, h = rhs_lukas(nu, rho, n, x, y, t)

    print('')
    print('           Minimum       Maximum')
    print('')
    print('  Ur:  %14.6g  %14.6g' % (np.min(f), np.max(f)))
    print('  Vr:  %14.6g  %14.6g' % (np.min(g), np.max(g)))
    print('  Pr:  %14.6g  %14.6g' % (np.min(h), np.max(h)))
#
#  Terminate.
#
    print('')
    print('rhs_lukas_test:')
    print('  Normal end of execution.')
    return


def rhs_poiseuille(nu, rho, n, x, y, t):

    f = np.zeros(n)
    g = np.zeros(n)
    h = np.zeros(n)

    return f, g, h


def rhs_poiseuille_test():

    nu = 1.0
    rho = 1.0

    print('')
    print('rhs_poiseuille_test')
    print('  Python version: %s' % (platform.python_version()))
    print('  Poiseuille Flow:')
    print('  Sample the Navier-Stokes right hand sides')
    print('  at the initial time T = 0, over a channel region.')
    print('  Kinematic viscosity NU = %g' % (nu))
    print('  Fluid density RHO = %g' % (rho))

    n = 1000
    x_lo = 0.0
    x_hi = 6.0
    y_lo = -1.0
    y_hi = +1.0
    x = x_lo + (x_hi - x_lo) * np.random.rand(n)
    y = y_lo + (y_hi - y_lo) * np.random.rand(n)
    t = 0.0

    f, g, h = rhs_poiseuille(nu, rho, n, x, y, t)

    print('')
    print('           Minimum       Maximum')
    print('')
    print('  Ur:  %14.6g  %14.6g' % (np.min(f), np.max(f)))
    print('  Vr:  %14.6g  %14.6g' % (np.min(g), np.max(g)))
    print('  Pr:  %14.6g  %14.6g' % (np.min(h), np.max(h)))
#
#  Terminate.
#
    print('')
    print('rhs_poiseuille_test:')
    print('  Normal end of execution.')
    return


def rhs_spiral(nu, rho, n, x, y, t):

    f = np.zeros(n)
    g = np.zeros(n)
    h = np.zeros(n)

    u = (1.0 + nu * t) * 2.0 \
        * (x ** 4 - 2.0 * x ** 3 + x ** 2) \
        * (2.0 * y ** 3 - 3.0 * y ** 2 + y)

    dudt = nu * 2.0 \
        * (x ** 4 - 2.0 * x ** 3 + x ** 2) \
        * (2.0 * y ** 3 - 3.0 * y ** 2 + y)

    dudx = (1.0 + nu * t) * 2.0 \
        * (4.0 * x ** 3 - 6.0 * x ** 2 + 2.0 * x) \
        * (2.0 * y ** 3 - 3.0 * y ** 2 + y)

    dudxx = (1.0 + nu * t) * 2.0 \
        * (12.0 * x ** 2 - 12.0 * x + 2.0) \
        * (2.0 * y ** 3 - 3.0 * y ** 2 + y)

    dudy = (1.0 + nu * t) * 2.0 \
        * (x ** 4 - 2.0 * x ** 3 + x ** 2) \
        * (6.0 * y ** 2 - 6.0 * y + 1.0)

    dudyy = (1.0 + nu * t) * 2.0 \
        * (x ** 4 - 2.0 * x ** 3 + x ** 2) \
        * (12.0 * y - 6.0)

    v = - (1.0 + nu * t) * 2.0 \
        * (2.0 * x ** 3 - 3.0 * x ** 2 + x) \
        * (y ** 4 - 2.0 * y ** 3 + y ** 2)

    dvdt = - nu * 2.0 \
        * (2.0 * x ** 3 - 3.0 * x ** 2 + x) \
        * (y ** 4 - 2.0 * y ** 3 + y ** 2)

    dvdx = - (1.0 + nu * t) * 2.0 \
        * (6.0 * x ** 2 - 6.0 * x + 1.0) \
        * (y ** 4 - 2.0 * y ** 3 + y ** 2)

    dvdxx = - (1.0 + nu * t) * 2.0 \
        * (12.0 * x - 6.0) \
        * (y ** 4 - 2.0 * y ** 3 + y ** 2)

    dvdy = - (1.0 + nu * t) * 2.0 \
        * (2.0 * x ** 3 - 3.0 * x ** 2 + x) \
        * (4.0 * y ** 3 - 6.0 * y ** 2 + 2.0 * y)

    dvdyy = - (1.0 + nu * t) * 2.0 \
        * (2.0 * x ** 3 - 3.0 * x ** 2 + x) \
        * (12.0 * y ** 2 - 12.0 * y + 2.0)

    p = rho * y
    dpdx = 0.0
    dpdy = rho

    f = dudt - nu * (dudxx + dudyy) \
        + u * dudx + v * dudy + dpdx / rho

    g = dvdt - nu * (dvdxx + dvdyy) \
        + u * dvdx + v * dvdy + dpdy / rho

    h = dudx + dvdy

    return f, g, h


def rhs_spiral_test():

    nu = 1.0
    rho = 1.0

    print('')
    print('rhs_spiral_test')
    print('  Python version: %s' % (platform.python_version()))
    print('  Spiral Flow:')
    print('  Sample the Navier-Stokes right hand sides')
    print('  at the initial time T = 0, over the unit square.')
    print('  Kinematic viscosity NU = %g' % (nu))
    print('  Fluid density RHO = %g' % (rho))

    n = 1000
    x_lo = 0.0
    x_hi = 1.0
    y_lo = 0.0
    y_hi = 1.0
    x = x_lo + (x_hi - x_lo) * np.random.rand(n)
    y = y_lo + (y_hi - y_lo) * np.random.rand(n)
    t = 0.0

    f, g, h = rhs_spiral(nu, rho, n, x, y, t)

    print('')
    print('           Minimum       Maximum')
    print('')
    print('  Ur:  %14.6g  %14.6g' % (np.min(f), np.max(f)))
    print('  Vr:  %14.6g  %14.6g' % (np.min(g), np.max(g)))
    print('  Pr:  %14.6g  %14.6g' % (np.min(h), np.max(h)))
#
#  Terminate.
#
    print('')
    print('rhs_spiral_test:')
    print('  Normal end of execution.')
    return


def rhs_taylor(nu, rho, n, x, y, t):

    f = np.zeros(n)
    g = np.zeros(n)
    h = np.zeros(n)

    return f, g, h


def rhs_taylor_test():

    nu = 1.0
    rho = 1.0

    print('')
    print('rhs_taylor_test')
    print('  Python version: %s' % (platform.python_version()))
    print('  Taylor Flow:')
    print('  Sample the Navier-Stokes right hand sides')
    print('  at the initial time T = 0, using a region that is')
    print('  the square centered at (1.5,1.5) with "radius" 1.0,')
    print('  Kinematic viscosity NU = %g' % (nu))
    print('  Fluid density RHO = %g' % (rho))

    n = 1000
    x_lo = 0.5
    x_hi = +2.5
    y_lo = 0.5
    y_hi = 2.5
    x = x_lo + (x_hi - x_lo) * np.random.rand(n)
    y = y_lo + (y_hi - y_lo) * np.random.rand(n)
    t = 0.0

    f, g, h = rhs_taylor(nu, rho, n, x, y, t)

    print('')
    print('           Minimum       Maximum')
    print('')
    print('  Ur:  %14.6g  %14.6g' % (np.min(f), np.max(f)))
    print('  Vr:  %14.6g  %14.6g' % (np.min(g), np.max(g)))
    print('  Pr:  %14.6g  %14.6g' % (np.min(h), np.max(h)))
#
#  Terminate.
#
    print('')
    print('rhs_taylor_test:')
    print('  Normal end of execution.')
    return


def rhs_vortex(nu, rho, n, x, y, t):

    f = - 2.0 * nu * (np.pi) ** 2 * (np.cos(np.pi * x) * np.sin(np.pi * y))
    g = 2.0 * nu * (np.pi) ** 2 * (np.sin(np.pi * x) * np.cos(np.pi * y))
    h = np.zeros(n)

    return f, g, h


def rhs_vortex_test():

    nu = 1.0
    rho = 1.0

    print('')
    print('rhs_vortex_test')
    print('  Python version: %s' % (platform.python_version()))
    print('  Sample the Vortex right hand sides')
    print('  at the initial time T = 0, using a region that is')
    print('  the square centered at (1.5,1.5) with "radius" 1.0,')
    print('  Kinematic viscosity NU = %g' % (nu))
    print('  Fluid density RHO = %g' % (rho))

    n = 1000
    x_lo = 0.5
    x_hi = +2.5
    y_lo = 0.5
    y_hi = 2.5
    x = x_lo + (x_hi - x_lo) * np.random.rand(n)
    y = y_lo + (y_hi - y_lo) * np.random.rand(n)
    t = 0.0

    f, g, h = rhs_vortex(nu, rho, n, x, y, t)

    print('')
    print('           Minimum       Maximum')
    print('')
    print('  Ur:  %14.6g  %14.6g' % (np.min(f), np.max(f)))
    print('  Vr:  %14.6g  %14.6g' % (np.min(g), np.max(g)))
    print('  Pr:  %14.6g  %14.6g' % (np.min(h), np.max(h)))
#
#  Terminate.
#
    print('')
    print('rhs_vortex_test:')
    print('  Normal end of execution.')
    return


def timestamp():

    t = time.time()
    print(time.ctime(t))

    return None


def uvp_gms(nu, rho, n, x, y, t):

    u, dudt, dudx, dudxx, dudy, dudyy, us, \
        v, dvdt, dvdx, dvdxx, dvdy, dvdyy, vs, \
        p, dpdt, dpdx, dpdxx, dpdy, dpdyy, ps = all_gms(nu, rho, n, x, y, t)

    return u, v, p


def uvp_gms_test():

    nu = 1.0
    rho = 1.0
    t = 1.0

    print('')
    print('uvp_gms_test')
    print('  Python version: %s' % (platform.python_version()))
    print('  GMS Flow:')
    print('  Estimate the range of velocity and pressure')
    print('  over the [-1,+1]x[-1,+1] square,')
    print('  at time T = 1.0')
    print('  Kinematic viscosity NU = %g' % (nu))
    print('  Fluid density RHO = %g' % (rho))

    n = 1000
    x_lo = -1.0
    x_hi = +1.0
    y_lo = -1.0
    y_hi = 1.0
    x = x_lo + (x_hi - x_lo) * np.random.rand(n)
    y = y_lo + (y_hi - y_lo) * np.random.rand(n)

    u, v, p = uvp_gms(nu, rho, n, x, y, t)

    print('')
    print('           Minimum       Maximum')
    print('')
    print('  U:  %14.6g  %14.6g' % (np.min(u), np.max(u)))
    print('  V:  %14.6g  %14.6g' % (np.min(v), np.max(v)))
    print('  P:  %14.6g  %14.6g' % (np.min(p), np.max(p)))
#
#  Terminate.
#
    print('')
    print('uvp_gms_test:')
    print('  Normal end of execution.')
    return


def uvp_gms_test2():

    r8_lo = -1.0
    r8_hi = +1.0

    nu = 1.0
    rho = 1.0
    t = 1.0

    print('')
    print('uvp_gms_test2')
    print('  Python version: %s' % (platform.python_version()))
    print('  GMS Flow:')
    print('  Estimate the range of velocity and pressure')
    print('  over the boundary of the [-1,+1]x[-1,+1] square,')
    print('  at time T = 1.0')
    print('  Kinematic viscosity NU = %g' % (nu))
    print('  Fluid density RHO = %g' % (rho))

    n = 400

    x = np.zeros(n)
    y = np.zeros(n)

    x[0:100] = np.linspace(r8_lo, r8_hi, 100)
    y[0:100] = r8_lo

    x[100:200] = r8_hi
    y[100:200] = np.linspace(r8_lo, r8_hi, 100)

    x[200:300] = np.linspace(r8_hi, r8_lo, 100)
    y[200:300] = r8_hi

    x[300:400] = r8_lo
    y[300:400] = np.linspace(r8_lo, r8_hi, 100)

    u, v, p = uvp_gms(nu, rho, n, x, y, t)

    print('')
    print('           Minimum       Maximum')
    print('')
    print('  U:  %14.6g  %14.6g' % (np.min(u), np.max(u)))
    print('  V:  %14.6g  %14.6g' % (np.min(v), np.max(v)))
    print('  P:  %14.6g  %14.6g' % (np.min(p), np.max(p)))
#
#  Terminate.
#
    print('')
    print('uvp_gms_test2:')
    print('  Normal end of execution.')
    return


def uvp_lukas(nu, rho, n, x, y, t):

    u = np.zeros(n)
    v = np.zeros(n)
    p = np.zeros(n)

    u = - np.cos(np.pi * x) / np.pi

    v = - y * np.sin(np.pi * x)

    return u, v, p


def uvp_lukas_test():

    nu = 1.0
    rho = 1.0

    print('')
    print('uvp_lukas_test')
    print('  Python version: %s' % (platform.python_version()))
    print('  Lukas Bystricky Flow:')
    print('  Estimate the range of velocity and pressure')
    print('  at the initial time T = 0, over the unit square.')
    print('  Kinematic viscosity NU = %g' % (nu))
    print('  Fluid density RHO = %g' % (rho))

    n = 1000
    x_lo = 0.0
    x_hi = +1.0
    y_lo = 0.0
    y_hi = 1.0
    x = x_lo + (x_hi - x_lo) * np.random.rand(n)
    y = y_lo + (y_hi - y_lo) * np.random.rand(n)
    t = 0.0

    u, v, p = uvp_lukas(nu, rho, n, x, y, t)

    print('')
    print('           Minimum       Maximum')
    print('')
    print('  U:  %14.6g  %14.6g' % (np.min(u), np.max(u)))
    print('  V:  %14.6g  %14.6g' % (np.min(v), np.max(v)))
    print('  P:  %14.6g  %14.6g' % (np.min(p), np.max(p)))
#
#  Terminate.
#
    print('')
    print('uvp_lukas_test:')
    print('  Normal end of execution.')
    return


def uvp_lukas_test2():

    r8_lo = 0.0
    r8_hi = +1.0

    nu = 1.0
    rho = 1.0
    t = 0.0

    print('')
    print('uvp_lukas_test2')
    print('  Python version: %s' % (platform.python_version()))
    print('  Lukas Bystricky Flow:')
    print('  Estimate the range of velocity and pressure')
    print('  on the boundary')
    print('  at the initial time T = 0, over the unit square.')
    print('  Kinematic viscosity NU = %g' % (nu))
    print('  Fluid density RHO = %g' % (rho))

    n = 400

    x = np.zeros(n)
    y = np.zeros(n)

    x[0:100] = np.linspace(r8_lo, r8_hi, 100)
    y[0:100] = r8_lo

    x[100:200] = r8_hi
    y[100:200] = np.linspace(r8_lo, r8_hi, 100)

    x[200:300] = np.linspace(r8_hi, r8_lo, 100)
    y[200:300] = r8_hi

    x[300:400] = r8_lo
    y[300:400] = np.linspace(r8_lo, r8_hi, 100)

    u, v, p = uvp_lukas(nu, rho, n, x, y, t)

    print('')
    print('           Minimum       Maximum')
    print('')
    print('  U:  %14.6g  %14.6g' % (np.min(u), np.max(u)))
    print('  V:  %14.6g  %14.6g' % (np.min(v), np.max(v)))
    print('  P:  %14.6g  %14.6g' % (np.min(p), np.max(p)))
#
#  Terminate.
#
    print('')
    print('uvp_lukas_test2:')
    print('  Normal end of execution.')
    return


def uvp_poiseuille(nu, rho, n, x, y, t):

    u = np.zeros(n)
    v = np.zeros(n)
    p = np.zeros(n)

    u = 1.0 - y ** 2
#
#  Can't write it this way or V becomes a scalar!
#
# v = 0.0;

    p = -2.0 * rho * nu * x

    return u, v, p


def uvp_poiseuille_test():

    nu = 1.0
    rho = 1.0

    print('')
    print('uvp_poiseuille_test')
    print('  Python version: %s' % (platform.python_version()))
    print('  Poiseuille Flow:')
    print('  Estimate the range of velocity and pressure')
    print('  at the initial time T = 0, over a channel region.')
    print('  Kinematic viscosity NU = %g' % (nu))
    print('  Fluid density RHO = %g' % (rho))

    n = 1000
    x_lo = +0.0
    x_hi = +6.0
    y_lo = -1.0
    y_hi = + 1.0
    x = x_lo + (x_hi - x_lo) * np.random.rand(n)
    y = y_lo + (y_hi - y_lo) * np.random.rand(n)
    t = 0.0

    u, v, p = uvp_poiseuille(nu, rho, n, x, y, t)

    print('')
    print('           Minimum       Maximum')
    print('')
    print('  U:  %14.6g  %14.6g' % (np.min(u), np.max(u)))
    print('  V:  %14.6g  %14.6g' % (np.min(v), np.max(v)))
    print('  P:  %14.6g  %14.6g' % (np.min(p), np.max(p)))
#
#  Terminate.
#
    print('')
    print('uvp_poiseuille_test:')
    print('  Normal end of execution.')
    return


def uvp_poiseuille_test2():

    x_lo = +0.0
    x_hi = +6.0
    y_lo = -1.0
    y_hi = + 1.0

    nu = 1.0
    rho = 1.0
    t = 0.0

    print('')
    print('uvp_poiseuille_test2')
    print('  Python version: %s' % (platform.python_version()))
    print('  Poiseuille Flow:')
    print('  Estimate the range of velocity and pressure')
    print('  on the boundary')
    print('  at the initial time T = 0, over a channel region.')
    print('  Kinematic viscosity NU = %g' % (nu))
    print('  Fluid density RHO = %g' % (rho))

    n = 400

    x = np.zeros(n)
    y = np.zeros(n)

    x[0:100] = np.linspace(x_lo, x_hi, 100)
    y[0:100] = y_lo

    x[100:200] = x_hi
    y[100:200] = np.linspace(y_lo, y_hi, 100)

    x[200:300] = np.linspace(x_hi, x_lo, 100)
    y[200:300] = y_hi

    x[300:400] = x_lo
    y[300:400] = np.linspace(y_hi, y_lo, 100)

    u, v, p = uvp_poiseuille(nu, rho, n, x, y, t)

    print('')
    print('           Minimum       Maximum')
    print('')
    print('  U:  %14.6g  %14.6g' % (np.min(u), np.max(u)))
    print('  V:  %14.6g  %14.6g' % (np.min(v), np.max(v)))
    print('  P:  %14.6g  %14.6g' % (np.min(p), np.max(p)))
#
#  Terminate.
#
    print('')
    print('uvp_poiseuille_test2:')
    print('  Normal end of execution.')
    return


def uvp_spiral(nu, rho, n, x, y, t):

    u = np.zeros(n)
    v = np.zeros(n)
    p = np.zeros(n)

    u = (1.0 + nu * t) * 2.0 \
        * x ** 2 * (x - 1.0) ** 2 \
        * y * (2.0 * y - 1.0) * (y - 1.0)

    v = - (1.0 + nu * t) * 2.0 \
        * x * (2.0 * x - 1.0) * (x - 1.0) \
        * y ** 2 * (y - 1.0) ** 2

    p = rho * y

    return u, v, p


def uvp_spiral_test():

    nu = 1.0
    rho = 1.0

    print('')
    print('uvp_spiral_test')
    print('  Python version: %s' % (platform.python_version()))
    print('  Spiral Flow:')
    print('  Estimate the range of velocity and pressure')
    print('  at the initial time T = 0, over the unit square.')
    print('  Kinematic viscosity NU = %g' % (nu))
    print('  Fluid density RHO = %g' % (rho))

    n = 1000
    x_lo = 0.0
    x_hi = +1.0
    y_lo = 0.0
    y_hi = 1.0
    x = x_lo + (x_hi - x_lo) * np.random.rand(n)
    y = y_lo + (y_hi - y_lo) * np.random.rand(n)
    t = 0.0

    u, v, p = uvp_spiral(nu, rho, n, x, y, t)

    print('')
    print('           Minimum       Maximum')
    print('')
    print('  U:  %14.6g  %14.6g' % (np.min(u), np.max(u)))
    print('  V:  %14.6g  %14.6g' % (np.min(v), np.max(v)))
    print('  P:  %14.6g  %14.6g' % (np.min(p), np.max(p)))
#
#  Terminate.
#
    print('')
    print('uvp_spiral_test:')
    print('  Normal end of execution.')
    return


def uvp_spiral_test2():

    r8_lo = 0.0
    r8_hi = +1.0

    nu = 1.0
    rho = 1.0
    t = 0.0

    print('')
    print('uvp_spiral_test2')
    print('  Python version: %s' % (platform.python_version()))
    print('  Spiral Flow:')
    print('  Estimate the range of velocity and pressure')
    print('  on the boundary')
    print('  at the initial time T = 0, over the unit square.')
    print('  Kinematic viscosity NU = %g' % (nu))
    print('  Fluid density RHO = %g' % (rho))

    n = 400

    x = np.zeros(n)
    y = np.zeros(n)

    x[0:100] = np.linspace(r8_lo, r8_hi, 100)
    y[0:100] = r8_lo

    x[100:200] = r8_hi
    y[100:200] = np.linspace(r8_lo, r8_hi, 100)

    x[200:300] = np.linspace(r8_hi, r8_lo, 100)
    y[200:300] = r8_hi

    x[300:400] = r8_lo
    y[300:400] = np.linspace(r8_lo, r8_hi, 100)

    u, v, p = uvp_spiral(nu, rho, n, x, y, t)

    print('')
    print('           Minimum       Maximum')
    print('')
    print('  U:  %14.6g  %14.6g' % (np.min(u), np.max(u)))
    print('  V:  %14.6g  %14.6g' % (np.min(v), np.max(v)))
    print('  P:  %14.6g  %14.6g' % (np.min(p), np.max(p)))
#
#  Terminate.
#
    print('')
    print('uvp_spiral_test2:')
    print('  Normal end of execution.')
    return


def uvp_taylor(nu, rho, n, x, y, t):

    cx = np.cos(np.pi * x)
    cy = np.cos(np.pi * y)
    c2x = np.cos(2.0 * np.pi * x)
    c2y = np.cos(2.0 * np.pi * y)
    sx = np.sin(np.pi * x)
    sy = np.sin(np.pi * y)
    e2t = np.exp(- 2.0 * np.pi * np.pi * nu * t)
    e4t = np.exp(- 4.0 * np.pi * np.pi * nu * t)

    u = np.zeros(n)
    v = np.zeros(n)
    p = np.zeros(n)

    u = - cx * sy * e2t
    v = sx * cy * e2t
    p = - 0.25 * rho * (c2x + c2y) * e4t

    return u, v, p


def uvp_taylor_test():

    nu = 1.0
    rho = 1.0

    print('')
    print('uvp_taylor_test')
    print('  Python version: %s' % (platform.python_version()))
    print('  Estimate the range of velocity and pressure')
    print('  at the initial time T = 0, using a region that is')
    print('  the square centered at (1.5,1.5) with "radius" 1.0,')
    print('  Kinematic viscosity NU = %g' % (nu))
    print('  Fluid density RHO = %g' % (rho))

    n = 1000
    x_lo = 0.5
    x_hi = +2.5
    y_lo = 0.5
    y_hi = 2.5
    x = x_lo + (x_hi - x_lo) * np.random.rand(n)
    y = y_lo + (y_hi - y_lo) * np.random.rand(n)
    t = 0.0

    u, v, p = uvp_taylor(nu, rho, n, x, y, t)

    print('')
    print('           Minimum       Maximum')
    print('')
    print('  U:  %14.6g  %14.6g' % (np.min(u), np.max(u)))
    print('  V:  %14.6g  %14.6g' % (np.min(v), np.max(v)))
    print('  P:  %14.6g  %14.6g' % (np.min(p), np.max(p)))
#
#  Terminate.
#
    print('')
    print('uvp_taylor_test:')
    print('  Normal end of execution.')
    return


def uvp_taylor_test2():

    r8_lo = 0.5
    r8_hi = +2.5

    nu = 1.0
    rho = 1.0
    t = 0.0

    print('')
    print('uvp_taylor_test2')
    print('  Python version: %s' % (platform.python_version()))
    print('  Estimate the range of velocity and pressure')
    print('  on the boundary')
    print('  at the initial time T = 0, using a region that is')
    print('  the square centered at (1.5,1.5) with "radius" 1.0,')
    print('  Kinematic viscosity NU = %g' % (nu))
    print('  Fluid density RHO = %g' % (rho))

    n = 400

    x = np.zeros(n)
    y = np.zeros(n)
#
#  Python is consistent in its willful flouting of sensible conventions.
#  X[0:100] means X from 0 to 99...!
#
    x[0:100] = np.linspace(r8_lo, r8_hi, 100)
    y[0:100] = r8_lo

    x[100:200] = r8_hi
    y[100:200] = np.linspace(r8_lo, r8_hi, 100)

    x[200:300] = np.linspace(r8_hi, r8_lo, 100)
    y[200:300] = r8_hi

    x[300:400] = r8_lo
    y[300:400] = np.linspace(r8_lo, r8_hi, 100)

    u, v, p = uvp_taylor(nu, rho, n, x, y, t)

    print('')
    print('           Minimum       Maximum')
    print('')
    print('  U:  %14.6g  %14.6g' % (np.min(u), np.max(u)))
    print('  V:  %14.6g  %14.6g' % (np.min(v), np.max(v)))
    print('  P:  %14.6g  %14.6g' % (np.min(p), np.max(p)))
#
#  Terminate.
#
    print('')
    print('uvp_taylor_test2:')
    print('  Normal end of execution.')
    return


def uvp_vortex(nu, rho, n, x, y, t):

    cx = np.cos(np.pi * x)
    cy = np.cos(np.pi * y)
    c2x = np.cos(2.0 * np.pi * x)
    c2y = np.cos(2.0 * np.pi * y)
    sx = np.sin(np.pi * x)
    sy = np.sin(np.pi * y)

    u = np.zeros(n)
    v = np.zeros(n)
    p = np.zeros(n)

    u = - cx * sy
    v = sx * cy
    p = - 0.25 * rho * (c2x + c2y)

    return u, v, p


def uvp_vortex_test():

    nu = 1.0
    rho = 1.0

    print('')
    print('uvp_vortex_test')
    print('  Python version: %s' % (platform.python_version()))
    print('  Sample the Vortex solution.')
    print('  Estimate the range of velocity and pressure')
    print('  at the initial time T = 0, using a region that is')
    print('  the square centered at (1.5,1.5) with "radius" 1.0,')
    print('  Kinematic viscosity NU = %g' % (nu))
    print('  Fluid density RHO = %g' % (rho))

    n = 1000
    x_lo = 0.5
    x_hi = +2.5
    y_lo = 0.5
    y_hi = 2.5
    x = x_lo + (x_hi - x_lo) * np.random.rand(n)
    y = y_lo + (y_hi - y_lo) * np.random.rand(n)
    t = 0.0

    u, v, p = uvp_vortex(nu, rho, n, x, y, t)

    print('')
    print('           Minimum       Maximum')
    print('')
    print('  U:  %14.6g  %14.6g' % (np.min(u), np.max(u)))
    print('  V:  %14.6g  %14.6g' % (np.min(v), np.max(v)))
    print('  P:  %14.6g  %14.6g' % (np.min(p), np.max(p)))
#
#  Terminate.
#
    print('')
    print('uvp_vortex_test:')
    print('  Normal end of execution.')
    return


def uvp_vortex_test2():

    r8_lo = 0.5
    r8_hi = +2.5

    nu = 1.0
    rho = 1.0
    t = 0.0

    print('')
    print('uvp_vortex_test2')
    print('  Python version: %s' % (platform.python_version()))
    print('  Sample the Vortex solution.')
    print('  Estimate the range of velocity and pressure')
    print('  on the boundary')
    print('  at the initial time T = 0, using a region that is')
    print('  the square centered at (1.5,1.5) with "radius" 1.0,')
    print('  Kinematic viscosity NU = %g' % (nu))
    print('  Fluid density RHO = %g' % (rho))

    n = 400

    x = np.zeros(n)
    y = np.zeros(n)
#
#  Python is consistent in its willful flouting of sensible conventions.
#  X[0:100] means X from 0 to 99...!
#
    x[0:100] = np.linspace(r8_lo, r8_hi, 100)
    y[0:100] = r8_lo

    x[100:200] = r8_hi
    y[100:200] = np.linspace(r8_lo, r8_hi, 100)

    x[200:300] = np.linspace(r8_hi, r8_lo, 100)
    y[200:300] = r8_hi

    x[300:400] = r8_lo
    y[300:400] = np.linspace(r8_lo, r8_hi, 100)

    u, v, p = uvp_vortex(nu, rho, n, x, y, t)
    return


if (__name__ == '__main__'):
    timestamp()
    navier_stokes_2d_exact_test()
    timestamp()
