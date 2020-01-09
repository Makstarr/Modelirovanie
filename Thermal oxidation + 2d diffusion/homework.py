from Oxidation import *
import Diff
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
Axes3D(fig)
surf, ax = plt.subplots()

ax = fig.gca(projection='3d')
Diff.zagonka()
surf = ax.plot_wireframe(Diff.Y1, Diff.X1, Diff.Z0, color="#727574")

temperature, orientation, p, xmax, B, A, ni = vvodDannih()

timeOxidation = 30*60  # cек
for t in range(1, timeOxidation):
    time += dt  # время которое обнуляется при изменении граничного элемента
    for j in range(0,n):
        # Окно
        if 2 < y[j] < 4:
            howMuchOxid[j] = X(time, A[j], B, xi[j]*2)/2
        # Левый заход за окно
        if y[j] <= 2:
            howMuchOxid[j] = (X(t, A[j], B, 0)/2*(1+erf(sqrt(2)/2*(y[j]-2) / (X(t, A[j], B, 0)))))
        # Правый заход за окно
        if y[j] >= 4:
            howMuchOxid[j] = (X(t, A[j], B, 0 * 2) / 2 * (1 + erf(sqrt(2) / 2 * (-(y[j] - 4)) / (X(t, A[j], B, 0 * 2)))))


    # Если самое большое значение оксида превзашло граничный элемент (плато)
    if max(howMuchOxid) >= borderElement:
        Diff.razgonka(time, max(howMuchOxid))  # Запускается разгонка от 0 до time (время окисления)
        borderElement += dx  # Граничный элемент следующий
        time -= time  # Время окисления обнуляется
        ax.plot(y, howMuchOxid, c="#e5d8c0")  # один из цветных графиков оксида (положительная часть)
        ax.plot(y, [-howMuchOxid[i] for i in range (0,n)], c="#e5d8c0")  # один из цветных графиков оксида (отрицательная часть)
        xi = copy.deepcopy(howMuchOxid)  # Копирует массив толщины для нового процесса
        temperature, orientation, p, xmax, B, A, ni = vvodDannih(Diff.Z0[:, j])

Diff.razgonka(time, max(howMuchOxid))  # Запускается разгонка от 0 до time (время окисления последнего кусочка)



surf = ax.plot_wireframe(Diff.Y1, Diff.X1, Diff.Z0, color="#141918")
oxid = ax.plot(y, howMuchOxid, c='black')
oxid2 = ax.plot(y, [-howMuchOxid[i] for i in range (0,n)], c='#5f4d3d')
ax.set_xlabel('Y' )
ax.set_ylabel('X')
ax.set_zlabel('Концентрация')
ax.set_ylim3d(-max(howMuchOxid)*1.2, 2)
Diff.pnFind(Diff.Z0)
ax.plot(Diff.x_pn,Diff.y_pn, c='red', label="pn-переход", linewidth=3)
plt.show()

