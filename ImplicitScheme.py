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


def b0(ht, hr):
    return 1 + 4 * gamma(ht, hr) + d * ht


def c0(ht, hr):
    return -4 * gamma(ht, hr)


def s0(u, ht):
    return u + ht * f_ri(0)


def b_i(hr, ht, ri):
    return 1 + 2 * gamma(ht, hr) + d * ht


def a_i(ht, hr, ri):
    return (gamma(ht, hr) * hr) / (2*ri) - gamma(ht, hr)


def c_i(ht, hr, ri):
    return -gamma(ht, hr) - (gamma(ht, hr) * hr)/(2*ri)


def s_i(u, ri, ht):
    return u + ht * f_ri(ri)

def a_I(ht, hr, ri):
    return -2*gamma(ht, hr)

def faster(args):
    hr, ht, riarr = args[0], args[1], args[2]
    res = [np.zeros(riarr.__len__())]
    len = riarr.__len__()
    coefficients_alpha_beta = [[0] * 2] * (len - 1)
    u = []
    values = np.zeros(len)

    for k in np.arange(1, int(180 / ht), 1):
        for i in np.arange(0, len, 1):
            if(i == 0):
                coefficients_alpha_beta[i] = [0, 0]
                values[i] = 0
            elif i == len - 1:
                values[i] = 0 # надо посчитать чему наверное
            else:
                if i == 1:
                    alpha0 = (-1 * c0(ht, hr)) / (b0(ht, hr))
                    beta0 = (s0(0,ht)) / (b0(ht, hr))
                    coefficients_alpha_beta[i] = [alpha0,beta0]

            elif i == len - 2:
            coefficients[i] = [-(((1 + 2 * gamma) * (1 + alp * h_x / k + c * h_x * h_x / (2 * k * h_t)) - gamma) / (
                        -gamma * coefficients[i - 1][0])),
                               ((values[i] + values[i] * (1 + 2 * gamma) * (c * h_x * h_x / (2 * k * h_t)) - (-gamma) *
                                 coefficients[i - 1][1]) / (-gamma * coefficients[i - 1][0]))]
        else:
            A = - gamma
            B = 1 + 2 * gamma
            C = - gamma
            F = values[i]
            Alpha = - (C/ (B + A * coefficients_alpha_beta[i - 1][0]))
            Betta = ((F - A * coefficients_alpha_beta[i - 1][1]) / (B + A * coefficients_alpha_beta[i - 1][0]))
           coefficients_alpha_beta[i] = [Alpha, Betta]
        for i in np.arange(len(x) - 2, -1, -1):
            values[i] = coefficients_alpha_beta[i][0] * values[i + 1] + coefficients_alpha_beta[i][1]

        return u



def xOy(args):
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
                ss.append([s_i(res[k - 1][j], ri, stept)])
            j += 1

        # правые части
        s = np.matrix(ss)

        # типо сетка
        u = np.zeros((len, len))

        # граничные условия
        u[0, 0] = b0(stept, stepr)
        u[0, 1] = c0(stept, stepr)
        u[len-1, len-1] = b_i(stepr, stept, riarr[len - 1])
        u[len-1, len-2] = a_I(stept, stepr, riarr[len - 1])

        # остальная система
        for i in np.arange(1, len-1, 1):
            u[i, i - 1] = a_i(stept, stepr, riarr[i])
            u[i, i] = b_i(stepr, stept, riarr[i])
            u[i, i + 1] = c_i(stept, stepr, riarr[i])

        temp = [y for y in np.linalg.solve(u, s).flat]
        res.append(temp)

    return res
