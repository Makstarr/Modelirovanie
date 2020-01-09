"""Разработать алгоритм и составить программу расчета ядерной и электронной
тормозных способностей ионов при ионной имплантации в соответствии с заданными условиями.
С помощью метода диффузионного приближения рассчитать проецированный пробег и разброс проецированных пробегов.
Для полученных значений рассчитать гауссовский профиль профиль распределения примеси."""

'''Провести расчеты ядерной и электронной тормозных способностей в соответствии с
заданными условиями. Рассчитать параметры Rp и Rp методом диффузионного приближения. 
Провести моделирования одномерного 
распределения примеси в заданном полупроводниковом 
материале после проведения процесса ионной имплантации.'''
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
qq = 14.4 # эВ*ангстрем
global Rp
Rp = 0
global dRpl

material = 'Si'
impurity = 'P'
epsilon = 11.7

E = 80 #Ev
Q = 1e15 #cm-3

Z2 = 14
Z1 = 15

M2 = 28
M1 = 31
N = 5e28  # m**-3 плотность атомов мишени

a0 = 0.529
ab = 0.8854*a0/((Z1**(2/3)+Z2**(2/3))**(1/2))
al = 0.8854*a0/((Z1**(2/3)+Z2**(2/3))**(1/2))
print (ab)
a = 1.1383
b = 0.01321
c = 0.21226
d = 0.19593

def Sne(E):
    e = ab*M2*E*1e3/(Z1*Z2*qq*(M1+M2))

    if e<10:
        return np.log(1+a*e)/2/(e+b**e**c+d*e**0.5)
    else:
        return np.log(e) / 2 / e

def Sn(E):
    return 8.462e-15*Z1*Z2*M1*Sne(E)/(M1+M2)/(Z1**0.23+Z2**0.23)*N*1e-6*1e-4*1e-3


def Se(E):
    k_LIH = Z1**(1/6)*0.0793*sqrt(Z1*Z2)*(M1+M2)**1.5/(Z1**(2/3)+Z2**(2/3))**0.75/(M1)**1.5/(M2)**0.5
    Cr = 2 * pi**1.5 * (ab * 1e-10) ** 2  * M1 * M2 / (M1 + M2) ** 2
    Ce =   ab * M2 / (Z1 * Z2 * qq/1e3 * (M1 + M2))
    K = k_LIH * Cr / (Ce**0.5)
    return (K*(E**0.5) * N )/1e6


def RdR(E):
    dE = 0.1
    Rp = dRpl = Rpsred =0
    ksi = 0
    for i in range(1, E*10):
        Rp = Rp*(1-M2/(2*M1)*Sn(i/10)/(Sn(i/10)+Se(i/10))*dE/(i/10))+dE/(Sn(i/10)+Se(i/10))
        ksi = ksi + 2*Rp/(Se(i/10)+Sn(i/10))*dE
        dRpl = dRpl + (ksi-2*dRpl)*M2/M1*Sn(i/10)/(Sn(i/10)+Se(i/10))*dE/(i/10)

    dRp = (ksi-Rp**2-dRpl)**0.5
    return Rp, dRp

def gauss(x, Rp, dRp):
    return Q/((2*pi)**0.5*dRp*1e-4)*exp(-(x - Rp)**2/(2*dRp**2))


def implant():
    sn = Sn(E)
    se = Se(E)
    print("Sn = %0.2f кэВ/мкм\n" % (Sn(E)),"Se = %0.2f кэВ/мкм\n" % (se))
    rp, drp = RdR(E)
    print("Rp = %.4f\n" % rp,"dRp = %.4f\n" % drp)

    n = 1000
    xmax = 0.5
    dx = xmax / n
    C = np.empty(n)
    x = np.empty(n)
    for i in range(0, n):
        if i == 0:
            x[i] = 0
        else:
            x[i] = x[i - 1] + dx
        C[i] = gauss(x[i],rp,drp)

    plt.plot(x, C, c='#6897bb', label='C(x)')
    plt.xlabel('x, мкм')
    plt.ylabel('Концентрация С')
    plt.legend()
    plt.xlim(0, x[-1])
    plt.show()

implant()
def graphsS():
    Sng = np.empty(400)
    Seg = np.empty(400)
    Eg = range(0,400)
    for E in Eg:
        if E == 0:
            Sng[E] = Seg[E] = 0
        Sng[E] = Sn(E)
        Seg[E] = Se(E)

    plt.plot(Eg,Sng, c='#6897bb', label='Sn(E)')
    plt.plot(Eg,Seg, c='#18573b', label='Se(E)')
    plt.xlabel('E, КэВ')
    plt.ylabel('Se/Sn, КэВ/мкм')
    plt.legend()
    plt.xlim(0, 400)
    plt.show()
graphsS()