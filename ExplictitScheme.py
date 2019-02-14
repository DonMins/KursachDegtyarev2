import math
import time
import numpy as np


# Система I(r)
def I_r(ri):
    import main as ex
    if(ri <= ex.parameter_array['a']):
        a = ex.parameter_array['P'] / (math.pi * ex.parameter_array['a'] ** 2)
        return a
    return 0

def f_ri(i) :
    import main as ex
    return (ex.parameter_array['betta'] * I_r(i)) / ex.parameter_array['c']

def psi(ht):
    import main as ex
    return (2 * ex.parameter_array['alf'] * ht) / (ex.parameter_array['c'] * ex.parameter_array['l'])


def gamma(ht, hr):
    import main as ex
    return (ex.parameter_array['k'] * ht) / (ex.parameter_array['c'] * hr ** 2)


def w0(w0, w1, ht, hr):
    p0 = 1 - 4*gamma(ht, hr) - psi(ht)
    return w0*p0 + w1*4*gamma(ht, hr) + f_ri(0)*ht


def wi(wm, wi, wp, ri, ht, hr):
    p0 = 1 - 2*gamma(ht, hr)- psi(ht)
    p1 = gamma(ht, hr) - gamma(ht, hr)*hr/(2*ri)
    p2 = gamma(ht, hr)*hr/(2*ri)+gamma(ht, hr)
    return wi*p0 + wm*p1 + wp*p2 + ht*f_ri(ri)

def wI(wI,wI_1,ht,hr,ri):
    p0 = 1 - 2 * gamma(ht, hr) - psi(ht)
    return wI * p0 + wI_1 * 2 * gamma(ht, hr) + f_ri(ri) * ht

def all(w, ht, hr, ris):
    len = ris.__len__()
    cur = [w0(w[0], w[1], ht, hr)]
    for i in np.arange(1, len-1, 1):
        cur.append(wi(w[i-1], w[i], w[i+1], ris[i], ht, hr))
    cur.append(wI(w[len-1],w[len-2],ht,hr,ris[len-1]))
    return cur

def xOy(args):

    stepr, stept, riarr = args[0], args[1], args[2]
    res = [np.zeros((riarr.__len__(), 1))]
    for k in np.arange(1, int(180 / stept), 1):
        res.append(all(res[k - 1], stept, stepr, riarr))

    return res