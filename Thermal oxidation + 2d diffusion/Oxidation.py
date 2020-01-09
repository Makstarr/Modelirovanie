import math, time, sys
import numpy as np
from numpy import sqrt, arctan, inf, real, log
from math import exp, pi
from scipy.special import erf
import copy
import Diff
# шаг по времени
dt = 1
# константы
m_dn = 1.08
m_dp = 0.59
q = 1.6e-19
k = 1.38e-23 / q  # Ev
C = 1e18  # Первый слой


time = 30*60#int(input("Введите время, мин: "))*60
temperature = 1050 + 273
orientation = 1
p = 2

# Массивы времени и толщины
t = np.empty(time)
x = np.empty(time)

# Расчет параметров кремния
Eg = 1.17 - 4.73e-4 / (temperature + 636) * temperature ** 2
Ei = Eg/2 - k*temperature/4
Nv = 4.82e15 * (m_dp ** 1.5) * (temperature ** 1.5)
Nc = 4.82e15 * (m_dn ** 1.5) * (temperature ** 1.5)
ni = ((Nc * Nv) ** 0.5) * math.exp(-1 * Eg / (2 * temperature * k))


# расчет коэффициентов
def Cminus(temperature):
    return exp((Ei+0.57-Eg)/k/temperature)

def Cplus(temperature):
    return exp((0.35-Ei)/k/temperature)

def Cminusminus(temperature):
    return exp((2*Ei+1.25-3*Eg)/k/temperature)

def Vl(temperature):
    return 2620*exp(-11/k/temperature)

def Vp(temperature):
    return 9.63e-16 * exp(2.83/k/temperature)

def Vn(temperature, n):
    return (1 + Cplus(temperature)*ni/n + Cminus(temperature)*n/ni + Cminusminus(temperature)*(n/ni)**2)/(1+Cminus(temperature)+Cplus(temperature)+Cminusminus(temperature))

# расчет констант при p = 1атм
def K(Ao, Ea, temperature):
    return Ao * exp(- Ea / k / temperature)

# толщина окисла x(t)
def X(t, A, B, xi):
    # постоянная времени
    # вместо xi при каждом запуске подставляется предыдущее значение x
    tau = (xi**2+A*xi)/B
    # толщина окисла x(t)
    x = (A/2)*(sqrt(1+(t+tau)/A**2*4*B)-1)
    return x
Co = [1e20 for i in range(0,50)]
def vvodDannih(Ci=copy.deepcopy(Co)):
    # ввод данных пользователем и перевод в СИ

    temperature = 1300#int(input("Введите температуру, C: ")) + 273
    p = 2#float(input("Введите давление, атм: "))
    xmax = 0.1#float(input("Введите толщину слоя, мкм: "))
    C = copy.deepcopy(Ci)#float(input("Введите концентрацию, см-3: "))

    orientation = 1.68

    # Расчет параметров кремния
    Eg = 1.17 - 4.73e-4 / (temperature + 636) * temperature ** 2
    Ei = Eg / 2 - k * temperature / 4
    Nv = 4.82e15 * (m_dp ** 1.5) * (temperature ** 1.5)
    Nc = 4.82e15 * (m_dn ** 1.5) * (temperature ** 1.5)
    ni = ((Nc * Nv) ** 0.5) * math.exp(-1 * Eg / (2 * temperature * k))

    T = temperature
    B = (K(7, 0.78, T) * (1 + Vp(T) * ni ** 0.22) * 2)/60
    # Линейная константа
    B_A = np.empty(50)
    A = np.empty(50)
    Csi = 6e16
    for i in range(0, 50):
        if C[i] <= Csi:
            B_A[i] = ((K(2.96e6, 2.05, T) * (p / orientation) * (1 + Vl(T) * (Vn(T, Csi) - 1)))) / 60
            A[i] = B / B_A[i]
        else:
            B_A[i] = ((K(2.96e6, 2.05, T) * (p / orientation) * (1 + Vl(T) * (Vn(T, C[i]) - 1)))) / 60
            A[i] = B / B_A[i]

    print(A)
    print('\n')
    return temperature, orientation, p, xmax, B, A, ni


n = 50
dx = 6/50
borderElement = dx
time = 0
howMuchOxid = [0 for i in range(0,n)]
xi = [0 for i in range(0,n)]
dy = 6/50
dt = 1
y = [dy*i for i in range(0,n)]

