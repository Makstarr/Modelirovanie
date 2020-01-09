import matplotlib.pyplot as plt
import math
import numpy as np
from scipy import special
'''
Программа для моделирования процессов диффузии примеси As в кремниевую подложку
Параметры задаются пользователем через командную строку
Расчитывается глубина залегания p-n перехода
ТарасовМД, БЭН-16-1
МИСиС 2019
'''
# установка параметров для графиков
# обычный масштаб
plt.figure(1)
plt.xlim(0, 1)
plt.title('График загонки и разгонки примеси в обычном масштабе')
plt.xlabel('Коорината, мкм')
plt.ylabel('Концентрация примеси')

# логарифмический масштаб
plt.figure(2)
plt.xlim(0, 1)
plt.ylim(0, 50)
plt.title('График загонки и разгонки примеси в логарифмическом масштабе')
plt.xlabel('Коорината, мкм')
plt.ylabel('Концентрация примеси')

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
n = 1000#int(input("Введите число точек: "))
Csi = 6e16#float(input("Введите концентрацию в подложке: "))

dx = 4/n
# пустые массивы
C = np.empty(n)
C_log = np.empty(n)
x = np.empty(n)
# граничные значения
x[0] = 0
C[0] = Co
C_log[0] = math.log(Co)


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

#
#
# Измененный вариант диффузии примеси
#

def zagonka():
    # массив координат
    for i in range(1, n):
        x[i] = x[i - 1] + dx
    # пустые массивы
    a = np.empty(n)
    d = np.empty(n)
    b = np.empty(n)
    r = np.empty(n)
    delta = np.empty(n)
    lyamda = np.empty(n)
    # граничные условия
    d[0] = 0
    a[0] = 1
    b[0] = 0
    r[0] = Co
    d[n - 1] = 0
    a[n - 1] = 1
    b[n - 1] = 0
    r[n - 1] = 0
    delta[0] = -d[0] / a[0]
    lyamda[0] = r[0] / a[0]

    # функция загонки при одном времени
    def c_time():
        # прямая прогонка
        for i in range(1, n):
            a[i] = -(2 + (dx ** 2 * 1e-8) / (Dif(0, temperature1, ni(temperature1)) * dt))
            r[i] = -(((dx ** 2 * 1e-8) * 0) / (Dif(0, temperature1, ni(temperature1)) * dt))
            b[i] = 1
            d[i] = 1
        for i in range(1, n):
            delta[i] = -d[i] / (a[i] + b[i] * delta[i - 1])
            lyamda[i] = (r[i] - b[i] * lyamda[i - 1]) / (a[i] + b[i] * delta[i - 1])
        C[-1] = lyamda[-1]
        # обратная прогонка
        for i in range(n - 2, -1, -1):
            C[i] = delta[i] * C[i + 1] + lyamda[i]

    # запуск функции загонки в интервале времени от 1 до времени разгонки (time1) с шагом dt
    for j in np.arange(1, time1, 1):
        c_time()

    # обычный масштаб
    plt.figure(1)
    plt.plot(x, C, label="Начальное распределение", c="black")

    # логарифмический масштаб
    plt.figure(2)
    plt.plot(x, C_log, label="Начальное распределение", c='black')


# задается функция разгонки
def razgonka():
    # пустые массивы
    a = np.empty(n)
    d = np.empty(n)
    b = np.empty(n)
    r = np.empty(n)
    delta = np.empty(n)
    lyamda = np.empty(n)
    # граничные условия
    d[0] = 1
    a[0] = -1
    b[0] = 0
    r[0] = 0
    d[n - 1] = 0
    a[n - 1] = 1
    b[n - 1] = 0
    r[n - 1] = 0
    delta[0] = -d[0] / a[0]
    lyamda[0] = r[0] / a[0]

    # функция разгонки при одном времени
    def c_time():
        # прямая прогонка
        for i in range(1, n):
            a[i] = -(2 + (dx ** 2 * 1e-8) / (Dif(C[i], temperature2,  ni(temperature2)) * dt))
            r[i] = -(((dx ** 2 * 1e-8) * C[i]) / (Dif(C[i], temperature2, ni(temperature2)) * dt))
            b[i] = 1
            d[i] = 1
        for i in range(1, n):
            delta[i] = -d[i] / (a[i] + b[i] * delta[i - 1])
            lyamda[i] = (r[i] - b[i] * lyamda[i - 1]) / (a[i] + b[i] * delta[i - 1])
        C[-1] = lyamda[-1]
        # обратная прогонка
        for i in range(n - 2, -1, -1):
            # переписывание массивов концентраций
            C[i] = delta[i] * C[i + 1] + lyamda[i]

    # запуск функции разгонки в интервале времени от 1 до времени разгонки (time2) с шагом dt
    for j in np.arange(1, time2, 1):
        c_time()

    # обычный масштаб
    plt.figure(1)
    plt.plot(x, C, label="После разгонки", c='green')

    # логарифмический масштаб
    plt.figure(2)
    plt.plot(x, C_log, label="После разгонки", c='green')


# запуск функций
zagonka()
razgonka()

# поиск p-n перехода
for i in range(1, n - 1):
    # Если концентрация стала меньше подложки выводится предыдущая координата
    if C[i] - Csi < 0:
        print("\n\nC примеси - C подложки ")
        print("в точке предыдущей точке: " + str(C[i-1] - Csi))
        print("в точке: " + str(C[i] - Csi))
        print("\nГлубина залегания PN перехода: "+str((x[i-1])))
        break

# показать графики
plt.show()
