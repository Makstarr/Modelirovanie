import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy import special
from matplotlib import cm


'''
Программа для моделирования процессов диффузии примеси As в кремниевую подложку
Параметры задаются пользователем через командную строку
Расчитывается глубина залегания p-n перехода
ТарасовМД, БЭН-16-1
МИСиС 2019
'''


# шаг по времени
dt = 1
# константы
m_dn = 1.08
m_dp = 0.59
q = 1.6e-19
kb = 1.38e-23
k = 1.38e-23 / q  # Ev
Co = 3e20

# ввод данных пользователем и перевод в СИ
print("\n\n\nПрограмма для моделирования процессов диффузии As в кремниевую подложку\n")
temperature1 = 1373#int(input("Введите температуру загонки, C: ")) + 273
temperature2 = 1273#int(input("Введите температуру разгонки, C: ")) + 273
time1 = 900#int(input("Введите время загонки, мин: "))*60
time2 = 1800#int(input("Введите время разгонки, мин: "))*60
#n = 100 int(input("Введите число точек: "))
Csi = 6e16#float(input("Введите концентрацию в подложке: "))

n = 30
xmax = 0.5
# пустые массивы
C = np.empty(n)
x = np.empty(n)
Y =np.empty(n)

fig = plt.figure()
Axes3D(fig)
ax = fig.gca(projection='3d')

# граничные значения
x[0] = 0
dx = xmax / n
width = 10
end = 20 + width

# заполнение массива концентраций
for i in range(n):
    if i == 0:
        x[i] = 0
        Y[i] = 0
    else:
        x[i] = x[i - 1] + dx
        Y[i] = Y[i - 1] + dx


# Массивы в матрицы
X1, Y1 = np.meshgrid(x, Y)
# Пустая матрица для концентраций
Z0 = X1*0
Z2 = X1*0
Z1 = X1*0
# функция для расчета коэффициента диффузии
def Dif(C, temperature, ni):
    return 0.66 * math.exp(-3.44 / (k * temperature)) + 12 * C / ni * math.exp(-4.05 / (k * temperature))


# функция для расчета собственной концентрации кремния
def ni(temperature):
    Eg = 1.17 - 4.73e-4 / (temperature + 636) * temperature ** 2
    Nv = 4.82e15 * (m_dp ** 1.5) * (temperature ** 1.5)
    Nc = 4.82e15 * (m_dn ** 1.5) * (temperature ** 1.5)
    ni = ((Nc * Nv) ** 0.5) * math.exp(-1 * Eg / (2 * temperature * k))
    return ni


# расчет ni для загонки и разгонки
ni2=ni(temperature2)
ni1=ni(temperature1)

#
#
# Диффузия
#

def zagonka_razgonka():
    # пустые массивы
    a = np.empty(n)
    d = np.empty(n)
    b = np.empty(n)
    r = np.empty(n)
    delta = np.empty(n)
    lyamda = np.empty(n)

    def c_time(C, C_1, C1, T, ni, condition):
        # граничные условия
        if condition == '1' or condition == '2' or condition == '3' or condition == '4':
            if condition == '1':
                r[0] = Co
                d[0] = 0
                a[0] = 1
            elif condition == '2':
                r[0] = 0
                d[0] = 0
                a[0] = 1
            elif condition == '3':
                d[0] = 1
                a[0] = -1
                r[0] = 0
            elif condition == '4':
                d[0] = 1
                a[0] = -1
                r[0] = 0
            b[0] = 0
            d[n - 1] = 0
            a[n - 1] = 1
            b[n - 1] = 0
            r[n - 1] = 0
            delta[0] = 0
            lyamda[0] = r[0] / a[0]
            if condition == '1' or condition == '3':
                for i in range(1, n, 1):
                    a[i] = -(2 + (dx ** 2 * 1e-8)/(Dif(C[i], T, ni) * dt))
                    r[i] = -(C1[i] - (2 - (dx ** 2 * 1e-8)/(Dif(C[i], T, ni)* dt))* C[i] + C_1[i])
                    b[i] = 1
                    d[i] = 1
                for i in range(0, n, 1):
                    delta[i] = -d[i] / (a[i] + b[i] * delta[i - 1])
                    lyamda[i] = (r[i] - b[i] * lyamda[i - 1]) / (a[i] + b[i] * delta[i - 1])
            if condition == '2'  or condition == '4':
                for i in range(0, n, 1):
                    a[i] = -(2 + (dx ** 2 * 1e-8)/(Dif(C[i], T, ni) * dt))
                    r[i] = -(C1[i] - (2 - (dx ** 2 * 1e-8)/(Dif(C[i], T, ni)* dt))* C[i] + C_1[i])
                    b[i] = 1
                    d[i] = 1
                for i in range(0, n, 1):
                    delta[i] = -d[i] / (a[i] + b[i] * delta[i - 1])
                    lyamda[i] = (r[i] - b[i] * lyamda[i - 1]) / (a[i] + b[i] * delta[i - 1])
            C[-1] = lyamda[-1]
            for i in range(n - 2, -1, -1):
                C[i] = delta[i] * C[i + 1] + lyamda[i]
            return C

    # Загонка
    for k in np.arange(1, time1, dt):
        print(k)
        for xx in range(1, n-1):
            if x[xx]<=0.3 and x[xx]>=0.1:
                Z0[xx] = c_time(Z0[xx],Z0[xx-1],Z0[xx+1], temperature1, ni1, '1')
            if x[xx]>=0.3 or x[xx]<=0.1:
                Z0[xx] = c_time(Z0[xx], Z0[xx - 1], Z0[xx + 1], temperature1, ni1, '2')
        for yy in range(1, n-1):
           Z0[:, yy] = c_time(Z0[:,yy], Z0[:,yy - 1], Z0[:,yy + 1], temperature1,  ni1, '2')
    surf = ax.plot_wireframe(Y1, X1, Z0, color="green")

    # Разгонка
    for k in np.arange(1, time2, dt):
        print (k)
        for xx in range(1, n-1):
            Z0[xx] = c_time(Z0[xx],Z0[xx-1],Z0[xx+1], temperature2, ni2, '3')
        for yy in range(1, n-1):
            Z0[:, yy] = c_time(Z0[:,yy], Z0[:,yy - 1], Z0[:,yy + 1], temperature2, ni2, '4')

    surf2 = ax.plot_wireframe(Y1, X1, Z0, color="red")
    print('\n')

zagonka_razgonka()
plt.show()
