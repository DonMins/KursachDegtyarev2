import matplotlib.pyplot as mpl
import time
import scipy.special as spc
import numpy as np
import scipy as sp
import plotly.offline as py

import ExplictitScheme
import ImplicitScheme

parameter_array = {
    'R': 4,
    'l': 0.5,
    'k': 0.011,
    'alf': 0.005,
    'c': 1.6,
    'betta': 0.008,
    'P': 40,
    'a': 0.8,
    'T': 180}

array_norm_N = {'0.1':7,
                '0.01':18,
                '0.001':52,
                '0.0001':137,
                '0.00001':347,
                '0.000001':877,
                '0.0000001':2202}

dataChanged=[0]
def updateParameret_array(dan):
    parameter_array['R'] = dan['R']
    parameter_array['l'] = dan['l']
    parameter_array['k'] = dan['k']
    parameter_array['alf'] = dan['alf']
    parameter_array['c'] = dan['c']
    parameter_array['betta'] = dan['betta']
    parameter_array['P'] = dan['P']
    parameter_array['a'] = dan['a']
    parameter_array['T'] = dan['T']

# надо как нибудь оптимизировать
bessel0 = [0, ]
temp = spc.jn_zeros(1, 100000)

for t in temp:
    bessel0.append(t)
##
def Bn(n):
    j0 = spc.j0(bessel0[n])
    j1 = spc.j1(bessel0[n] / 5)
    numerator = 2 * parameter_array['betta'] * parameter_array['P'] * j1
    denominator = 5 * parameter_array['c'] * np.pi * (parameter_array['a']) ** 2 * bessel0[n] * (j0) ** 2
    return numerator/denominator

def An(n,t):
    numerator = Bn(n) * parameter_array['c'] * parameter_array['l'] * (parameter_array['R']) ** 2
    denominator = (2 * parameter_array['alf'] * (parameter_array['R']) ** 2) + (parameter_array['k'] * \
                                                                                bessel0[n] ** 2 * parameter_array['l'])

    x = (((-2 * parameter_array['alf'] * parameter_array['R'] ** 2) - (parameter_array['k'] * \
                                                                       bessel0[n] ** 2 * parameter_array['l'])) * t) / (parameter_array['c'] * parameter_array['l'] * parameter_array['R'] ** 2)

    all = numerator/denominator
    return all * (1- np.exp(x))

def analitic_eps(rank):
    N = ((2 * parameter_array['betta'] * parameter_array['R'] * parameter_array['P'] * 10 ** (1 / 2)) / (15 * (parameter_array['a']) ** 2 * parameter_array['k'] *
                                                                                                         (np.pi) ** 3 * rank))**(2/3)
    return round(N)

flag=[]
NMAX=[0]


def experimental_eps(r, t, eps):
    numerator = (parameter_array['l'] * parameter_array['betta'] * parameter_array['P'])
    denominator = (50 * sp.pi * parameter_array['alf'] * parameter_array['a'] ** 2)
    all = numerator / denominator
    x = ((-2 * parameter_array['alf'] * t) / (parameter_array['c'] * parameter_array['l']))
    u = all * (1 - np.exp(x))

    N = analitic_eps(eps)

    for i in np.arange(1, N+1, 1):
        u += An(i, t) * spc.j0((bessel0[i] * r) / parameter_array['R'])
        flag.append(u)

    i=N-2
    analit_eps = flag[N-1]
    while (abs ( analit_eps - flag[i])<=eps)and(i>0):
        i=i-1
    if ((i+2)>NMAX[0]):
        NMAX[0] = i+2
    flag.clear()
    return NMAX[0]

def u(r,t, eps,colslag):

    numerator = (parameter_array['l'] * parameter_array['betta'] * parameter_array['P'])
    denominator = (50 * sp.pi * parameter_array['alf'] * parameter_array['a'] ** 2)
    all = numerator/denominator
    x = ((-2 * parameter_array['alf'] * t) / (parameter_array['c'] * parameter_array['l']))
    u = all * (1- np.exp(x))
    if (colslag==0):
        if dataChanged[0]==1:

            N = analitic_eps(eps)
            for i in np.arange(1, N, 1):
                u += An(i, t) * spc.j0((bessel0[i] * r) / parameter_array['R'])

        else:
            N = array_norm_N[str(eps)]

            for i in np.arange(1,N,1):
                u +=An(i, t) * spc.j0((bessel0[i] * r) / parameter_array['R'])

    else:
        for i in np.arange(1,colslag,1):
            u +=An(i, t) * spc.j0((bessel0[i] * r) / parameter_array['R'])

    return u

def getparams():
    return parameter_array



if __name__ == '__main__':
    exit = 0
    while(exit==0):
        parameter_array = getparams()
        print("Текущие параметры")
        print("R = {0}  l = {1}  k = {2}".format(parameter_array['R'], parameter_array['l'], parameter_array['k']))
        print("alfa = {0} c = {1} betta = {2}".format(parameter_array['alf'], parameter_array['c'], parameter_array['betta']))
        print("P = {0}, a = {1}".format(parameter_array['P'], parameter_array['a']))
        yes = input("Введите +, если хотите изменить параметры  иначе - ")
        if yes.upper() == "+":
            dataChanged[0]=1
            while True:
                try:
                    R = float(input("Введите R "))
                    if((R<=0) or R>=10):
                        print("Вы должны ввести положительное число больше, 0 и меньше 10, попробуйте снова.")
                    else:
                        break
                except ValueError:
                    print("Вы должны ввести положительное число больше, 0 и меньше 10, попробуйте снова.")

            while True:
                try:
                    l = float(input("Введите l "))
                    if  (l<=0 or l>=5):
                        print("Вы должны ввести положительное число, больше 0 и меньше 5 , попробуйте снова.")
                    else:
                        break
                except ValueError:
                    print("Вы должны ввести положительное число, больше 0 и меньше 5 , попробуйте снова.")

            while True:
                try:
                    k = float(input("Введите k "))
                    if  ((k<=0) or (k>1)):
                        print("Вы должны ввести положительное число, больше 0 и меньше 1 , попробуйте снова.")
                    else:
                        break
                except ValueError:
                    print("Вы должны ввести положительное число, больше 0 и меньше 1 , попробуйте снова.")

            while True:
                try:
                    alfa = float(input("Введите alfa "))
                    if((alfa<=0) or (alfa>1)):
                        print("Вы должны ввести положительное число, больше 0 и меньше 1 , попробуйте снова.")
                    else:
                        break
                except ValueError:
                    print("Вы должны ввести положительное число, больше 0 и меньше 1 , попробуйте снова.")

            while True:
                try:
                    c = float(input("Введите c "))
                    if ((c<=0) or (c>10)):
                        print("Вы должны ввести положительное число, больше 0 и меньше 10 , попробуйте снова.")
                    else:
                        break
                except ValueError:
                    print("Вы должны ввести положительное число, больше 0 и меньше 10 , попробуйте снова.")

            while True:
                try:
                    betta = float(input("Введите betta "))
                    if ((betta<=0) or (betta>1)):
                        print("Вы должны ввести положительное число, больше 0 и меньше 1 , попробуйте снова.")
                    else:
                        break
                except ValueError:
                    print("Вы должны ввести положительное число, больше 0 и меньше 1 , попробуйте снова.")

            while True:
                try:
                    P = float(input("Введите P "))
                    if(P<=0):
                        print("Вы должны ввести положительное число, попробуйте снова.")
                    else:
                        break
                except ValueError:
                    print("Вы должны ввести положительное число, попробуйте снова.")

            while True:
                try:
                    a = float(input("Введите a "))
                    if(a <= 0):
                        print("Вы должны ввести целое положительное число , попробуйте снова.")
                    else:
                        break
                except ValueError:
                    print("Вы должны ввести целое положительное число , попробуйте снова.")


            parameter_array['R'] = float(R)
            parameter_array['l'] = float(l)
            parameter_array['k'] = float(k)
            parameter_array['alf'] = float(alfa)
            parameter_array['c'] = float(c)
            parameter_array['betta'] = float(betta)
            parameter_array['P'] = float(P)
            parameter_array['a'] = float(a)

            updateParameret_array(parameter_array)
            parameter_array = getparams()
            print("Текущие параметры")
            print("R = {0} l = {1} k = {2}".format(parameter_array['R'], parameter_array['l'], parameter_array['k']))
            print("alfa = {0} c = {1} betta = {2}".format(parameter_array['alf'], parameter_array['c'], parameter_array['betta']))
            print("P = {0}, a = {1}".format(parameter_array['P'], parameter_array['a']))

        while True:
            try:
                a = int(input("Введите количество шагов по R "))
                if(a <= 0):
                    print("Вы должны ввести целое положительное число , попробуйте снова.")
                else:
                    break
            except ValueError:
                print("Вы должны ввести целое положительное число , попробуйте снова.")

        stepr = parameter_array['R']/a

        while True:
            try:
                a = int(input('и по T '))
                if(a <= 0):
                    print("Вы должны ввести целое положительное число , попробуйте снова.")
                else:
                    break
            except ValueError:
                print("Вы должны ввести целое положительное число , попробуйте снова.")

        stept = parameter_array['T']/a
        while True:
            try:
                a = int(input("Введите момент времени "))
                if(a <= 0):
                    print("Вы должны ввести целое положительное число , попробуйте снова.")
                else:
                    break
            except ValueError:
                print("Вы должны ввести целое положительное число , попробуйте снова.")
        curtime = a
        print("Подождите, пожалуйста, выполняются вычисления...")

        riarr = [st for st in np.arange(0, parameter_array['R'], stepr)]
        args =(stepr, stept, riarr)
        t1= time.time()

        impicit = ImplicitScheme.xOy(args)

        print('Время работы неявной схемы' + ' {0:.2f}'.format(time.time() - t1))

        t2 = time.time()

            #explicit = ExplictitScheme.xOy(args)
        print('Время работы явной схемы' + ' {0:.2f}'.format(time.time() - t2))


        y1 = [u(step,curtime, 0.01,0) for step in riarr]
        y2 = impicit[int(curtime / stept)]
        #y2 = explicit[int(curtime / stept)]
        ln0, ln1 = mpl.plot(riarr,y2,riarr,y1)
        mpl.legend((ln0, ln1), ('явная', 'Аналитическое', "Явная"),
                   title='R: {0}, l: {1}, k: {2}, alf: {3}, c: {4}, betta: {5}, P: {6}, a: {7} \n step for R : {8} \n step '
                         'for T: {9} \n time = {10}'.format(parameter_array['R'], parameter_array['l'], parameter_array['k'],
                                                            parameter_array['alf'], parameter_array['c'],
                                                            parameter_array['betta'],
                                                            parameter_array['P'], parameter_array['a'], stepr, stept, curtime))

        mpl.xlabel('Радиус')
        mpl.ylabel('Tемпература')
        mpl.grid()
        mpl.show()
        exit = int(input("Для выхода нажмите - 1 , чтобы продолжить - 0"))
        dataChanged[0]=0
    print("-------------------------------------------------------------------------------------")
