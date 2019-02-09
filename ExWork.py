import scipy.special as spc
import numpy as np
import scipy as sp
import plotly.offline as py
import numbers
import matplotlib.pyplot as plt




# Для построения графика погрешности
Xnew = []
Ynew=[]

def plot_eps(x,y):

    fig=  plt.figure("График погрешности ")
    ax = fig.add_subplot(111)
    plt.grid(True)
    ax.set_yticks([2202,1900,1700,1500,1250,1000,877,650,500,347,137,7])
    plt.plot(x, y, marker='o',markersize = 2.5)
    leg1, leg2= plt.plot(x, y,Xnew,Ynew,linewidth = 1)
    plt.legend((leg1, leg2), ("График достаточного числа членов ряда", "Интерполированный график достаточного числа членов ряда "))
    #plt.plot(Xnew,Ynew,color='#FF0000')
    plt.xlabel("accuracy($\epsilon$)")
    plt.ylabel("number of members (N)")

    ax.set_xscale('log')



    plt.show()

#Аппроксимация графика погрешности
# X = [0.0000001,0.000001,0.00001,0.0001,0.001,0.01,0.1]
# Y=[2202,877,347,137,52,18,7]
#
# yt = 0.0000001
# st = -8
# i =0
# while (yt<=10**-1):
#     Ynew.append(experimental_eps(0,180,yt,parameter_array))
#     NMAX[0]=0
#     Xnew.append(yt)
#     yt= yt + 9*10**st
#     i=i+1
#     if(i==9):
#         st=st+1
#         i=0
#         yt=10**(st+1)
#
#
# plot_eps(X,Y)
#
#
# exit = 0
# while(exit==0):
#     parameter_array = getparams()
#     print("Текущие параметры")
#     print("R = {0} l = {1} k = {2}".format(parameter_array['R'], parameter_array['l'], parameter_array['k']))
#     print("alfa = {0} c = {1} betta = {2}".format(parameter_array['alf'], parameter_array['c'], parameter_array['betta']))
#     print("P = {0}, a = {1}, T={2}".format(parameter_array['P'], parameter_array['a'], parameter_array['T']))
#     yes = input("Введите +, если хотите изменить параметры  иначе - ")
#
#     if yes.upper() == "+":
#         dataChanged[0]=1
#         while True:
#             try:
#                 R = float(input("Введите R "))
#                 if((R<=0) or R>=10):
#                     print("Вы должны ввести положительное число больше, 0 и меньше 10, попробуйте снова.")
#                 else:
#                     break
#             except ValueError:
#                 print("Вы должны ввести положительное число больше, 0 и меньше 10, попробуйте снова.")
#
#         while True:
#             try:
#                 l = float(input("Введите l "))
#                 if  (l<=0 or l>=5):
#                     print("Вы должны ввести положительное число, больше 0 и меньше 5 , попробуйте снова.")
#                 else:
#                      break
#             except ValueError:
#                 print("Вы должны ввести положительное число, больше 0 и меньше 5 , попробуйте снова.")
#
#         while True:
#             try:
#                 k = float(input("Введите k "))
#                 if  ((k<=0) or (k>1)):
#                     print("Вы должны ввести положительное число, больше 0 и меньше 1 , попробуйте снова.")
#                 else:
#                     break
#             except ValueError:
#                 print("Вы должны ввести положительное число, больше 0 и меньше 1 , попробуйте снова.")
#
#         while True:
#             try:
#                 alfa = float(input("Введите alfa "))
#                 if((alfa<=0) or (alfa>1)):
#                     print("Вы должны ввести положительное число, больше 0 и меньше 1 , попробуйте снова.")
#                 else:
#                     break
#             except ValueError:
#                 print("Вы должны ввести положительное число, больше 0 и меньше 1 , попробуйте снова.")
#
#         while True:
#             try:
#                 c = float(input("Введите c "))
#                 if ((c<=0) or (c>10)):
#                     print("Вы должны ввести положительное число, больше 0 и меньше 10 , попробуйте снова.")
#                 else:
#                     break
#             except ValueError:
#                 print("Вы должны ввести положительное число, больше 0 и меньше 10 , попробуйте снова.")
#
#         while True:
#             try:
#                 betta = float(input("Введите betta "))
#                 if ((betta<=0) or (betta>1)):
#                     print("Вы должны ввести положительное число, больше 0 и меньше 1 , попробуйте снова.")
#                 else:
#                     break
#             except ValueError:
#                 print("Вы должны ввести положительное число, больше 0 и меньше 1 , попробуйте снова.")
#
#         while True:
#             try:
#                 P = float(input("Введите P "))
#                 if(P<=0):
#                     print("Вы должны ввести положительное число, попробуйте снова.")
#                 else:
#                     break
#             except ValueError:
#                 print("Вы должны ввести положительное число, попробуйте снова.")
#
#         while True:
#             try:
#                 T = int(input("Введите T "))
#                 if(T <= 0):
#                      print("Вы должны ввести целое положительное число , попробуйте снова.")
#                 else:
#                     break
#             except ValueError:
#                 print("Вы должны ввести целое положительное число , попробуйте снова.")
#
#         while True:
#             try:
#                 a = float(input("Введите a "))
#                 if(a <= 0):
#                     print("Вы должны ввести целое положительное число , попробуйте снова.")
#                 else:
#                      break
#             except ValueError:
#                 print("Вы должны ввести целое положительное число , попробуйте снова.")
#
#
#         parameter_array['R'] = float(R)
#         parameter_array['l'] = float(l)
#         parameter_array['k'] = float(k)
#         parameter_array['alf'] = float(alfa)
#         parameter_array['c'] = float(c)
#         parameter_array['betta'] = float(betta)
#         parameter_array['P'] = float(P)
#         parameter_array['a'] = float(a)
#         parameter_array['T'] = float(T)
#         updateParameret_array(parameter_array)
#         parameter_array = getparams()
#         print("Текущие параметры")
#         print("R = {0} l = {1} k = {2}".format(parameter_array['R'], parameter_array['l'], parameter_array['k']))
#         print("alfa = {0} c = {1} betta = {2}".format(parameter_array['alf'], parameter_array['c'], parameter_array['betta']))
#         print("P = {0}, a = {1}, T={2}".format(parameter_array['P'], parameter_array['a'], parameter_array['T']))
#
#     scoreOrSeries = (input("Хотите вычислить ряд ( введите - 1) или погрешность (введите 2) ? "))
#
#     if scoreOrSeries == "1":
#         yes_no = (input(
#             "Хотите посчитать ряд используя точность( введите - 1) или количество слагаемых для посчета(введите 2) ? "))
#
#         if yes_no == "1":
#
#             while True:
#                 try:
#                     eps = float(input("Введите точность (пример 0.001 или 10e-4)"))
#                     if ((eps<0) or (eps>=0.9)):
#                         print("Вы должны ввести число больше 0 и меньше 0.9 , попробуйте снова.")
#                     else:
#                         break
#                 except ValueError:
#                     print("Вы должны ввести число больше 0 и меньше 0.9 , попробуйте снова.")
#
#             print("Подождите, пожалуйста, выполняются вычисления...")
#             drawUIR(parameter_array, eps, 0)
#             drawUIT(parameter_array, eps, 0)
#
#         if yes_no == "2":
#
#             while True:
#                 try:
#                     n = int(input("Введите количество слагаемых N = "))
#                     if ((n<=3) or (n>10000)):
#                         print("Вы должны ввести целое число > 3 и < 10 000, попробуйте снова.")
#                     else:
#                         break
#                 except ValueError:
#                     print("Вы должны ввести целое число > 3 и < 10 000, попробуйте снова.")
#
#             print("Подождите, пожалуйста, выполняются вычисления...")
#             drawUIR(parameter_array, 0, n)
#             drawUIT(parameter_array, 0, n)
#
#     if scoreOrSeries == "2":
#         error = (input(
#             "Рассчитать кол -во слагаемых для достижени конкретной точности (введите 1 ) или для всего ряда (введите 2 (это займет много времени ))"))
#         if error == "1":
#             while True:
#                 try:
#                     ep = float(input("Введите точность (пример 0.001 или 10e-4)"))
#                     if ((ep<0) or (ep>=0.9)):
#                         print("Вы должны ввести число больше 0 и меньше 0.9 , попробуйте снова.")
#                     else:
#                         break
#                 except ValueError:
#                     print("Вы должны ввести число больше 0 и меньше 0.9 , попробуйте снова.")
#
#             print("Подождите, пожалуйста, выполняются вычисления...")
#             for t in np.arange(1, parameter_array['T'] + 1):
#                 for r in np.arange(0, parameter_array['R'] + 0.1, 0.1):
#                     experimental_eps(r, t, ep, parameter_array)
#             print("NMAX избыточое для eps = 10^-", ep, "равно ", analitic_eps(ep))
#             print("NMAX достаточное для eps = 10^-", ep, "равно ", NMAX)
#             NMAX[0] = 0
#             flag.clear()
#
#         if error == "2":
#             print("Подождите, пожалуйста, выполняются вычисления...")
#             Y = []
#             EPS = [0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1] # степени
#             for ep in EPS:
#                 for t in np.arange(1, parameter_array['T'] + 1):
#                     for r in np.arange(0, parameter_array['R'] + 0.1, 0.1):
#                         experimental_eps(r, t, ep, parameter_array)
#                 print("N избыточое для eps = ", ep, "равно ", analitic_eps(ep))
#                 print("N достаточное для eps = ", ep, "равно ", NMAX)
#                 Y.append(NMAX[0])
#                 NMAX[0] = 0
#                 flag.clear()
#
#             plot_eps(EPS, Y)
#     exit = int(input("Для выхода нажмите - 1 , чтобы продолжить - 0"))
#     dataChanged[0]=0
#     print("-------------------------------------------------------------------------------------")