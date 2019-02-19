import math
import time
import numpy as np

import main
parameter_array = main.getparams()
# Система I(r)
def I_r(ri):

    if (ri <= parameter_array['a']):
        a = parameter_array['P'] / (math.pi * parameter_array['a'] ** 2)
        return a
    return 0

# Элементы из первой замены
d = (2 * parameter_array['alf']) / (parameter_array['c'] * parameter_array['l'])

g = parameter_array['k'] / parameter_array['c']


def f_ri(i):
    return (parameter_array['betta'] * I_r(i)) / parameter_array['c']


def gamma(ht, hr):
    return g * ht / (hr ** 2)


def p0(ht, hr):
    return 1 + 4 * gamma(ht, hr) + d * ht


def q0(ht, hr):
    return 4 * gamma(ht, hr)


def s0(u, ht):
    return u + ht * f_ri(0)


def pi(hr, ht, ri):
    return 1 + 2 * gamma(ht, hr) + d * ht


def hi(ht, hr, ri):
    return gamma(ht, hr) - (gamma(ht, hr) * hr) / (2*ri)


def qi(ht, hr, ri):
    return gamma(ht, hr) + (gamma(ht, hr) * hr )/( 2*ri)


def si(u, ri, ht):
    return u + ht * f_ri(ri)

def hI(ht, hr, ri):
    return 2*gamma(ht, hr)


def xOy(args):
    yu = 0

    stepr, stept, riarr = args[0], args[1], args[2]
    res = [np.zeros(riarr.__len__())]
    len = riarr.__len__()
    for k in np.arange(1, int(180 / stept), 1):

        ss = []
        j = 0
        for ri in riarr[0:len]:
            if (ri == 0):
                ss.append([s0(0, stept)])
            else:
                ss.append([si(res[k - 1][j], ri, stept)])
            j += 1

        # правые части
        s = np.matrix(ss)

        # типо сетка
        u = np.zeros((len, len))

        # граничные условия
        u[0, 0] = p0(stept, stepr)
        u[0, 1] = -1 * q0(stept, stepr)
        u[len-1, len-1] = pi(stepr, stept, riarr[len-1])
        u[len-1, len-2] = -1*hI(stept, stepr, riarr[len-1])

        # остальная система
        for i in np.arange(1, len-1, 1):
            u[i, i - 1] = -1*hi(stept, stepr, riarr[i])
            u[i, i] = pi(stepr, stept, riarr[i])
            u[i, i + 1] = -1*qi(stept, stepr, riarr[i])

        if (yu==0):
            print(u)

        temp = [y for y in np.linalg.solve(u, s).flat]
        res.append(temp)
        temp=[]
        del s
        if (yu == 0):
            print(res)
            yu=1



    return res
