import matplotlib.pyplot as plt
import math, time, sys
import numpy as np
from numpy import sqrt, arctan, inf, real, log
from math import exp, pi

# шаг по времени
dt = 1
# константы
m_dn = 1.08
m_dp = 0.59
q = 1.6e-19
k = 1.38e-23 / q  # Ev
Co = 1e20

# ввод данных пользователем и перевод в СИ
print("\n\n\nПрограмма для моделирования окисления кремния\n")
temperature = int(input("Введите температуру, C: ")) + 273
time = int(input("Введите время, мин: "))*60
Csi = float(input("Введите концентрацию в подложке: "))
orientation = int(input("Введите ориентацию, <...> : "))
p = float(input("Введите давление, атм: "))

# Учет зависимости линейной константы от ориентации
if orientation == 100:
    orientation = 1.68
else:
    orientation = 1

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
def X(t, A, B, xi=0):
    # постоянная времени
    tau = (xi**2+A*xi)/B
    # толщина окисла x(t)
    x = (A/2)*(sqrt(1+(t+tau)/A**2*4*B)-1)
    return x

def okislenie(T, time, Csi):
    # Параболическая константа
    B = K(7, 0.78, T)*(1 + Vp(T)*ni**0.22)*2
    print('B = %.3f' % B)
    # Линейная константа
    B_A = (K(2.96e6, 2.05, T)*(p/orientation)*(1 + Vl(T)*(Vn(T, Csi)-1)))
    print('B/A = %.3f' % B_A)
    # Выражение А
    A=B/B_A
    print('A = %.3f' % A)

    # Массивы времени и толщины
    t = np.empty(time)
    x = np.empty(time)

    # Цикл расчета х(t)
    for i in range(0, time):
        if i == 0:
            t[i] = 0
        else:
            t[i] = t[i-1] + dt
        x[i] = X(t[i], A, B, 0)

    # Параметры графика
    plt.plot(t, x, c='black')
    plt.ylabel('Толщина слоя x, мкм')
    plt.xlabel('Время t, c')
    plt.xlim(0,time)
    plt.ylim(0,x[-1])
    plt.show()

okislenie(temperature, time, Csi)
