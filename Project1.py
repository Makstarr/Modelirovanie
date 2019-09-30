import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
fig = plt.figure()
fig.set_dpi(100)
ax1 = fig.add_subplot(1,1,1)

def f(x):
    return 1e+19*np.exp(-4*x)


def Dif(x):
    Dif =0.000000000001
    return Dif
xmax=4
n=300


a=np.empty(n+1)
b=np.empty(n+1)
d=np.empty(n+1)
r=np.empty(n+1)
delta=np.empty(n+1)
lambd=np.empty(n+1)

dx = xmax/(n)
x=np.empty(n+1)
c=np.empty(n+1)

c[0] = 1e+19
x[0] = 0

x[-1]=xmax
c[n]=0
d[n]=0
a[n]=1
b[n]=0
r[n]=0
d[1]=1
a[1]=-1
b[1]=0
r[1]=0
delta[1]=-d[1]/a[1]
lambd[1]=r[1]/a[1]
cc=[]
dt=1

if __name__ == '__main__':
    for i in range(0, n ):
        x[i+1] = x[i] + dx
        c[i+1] = f(x[i+1])
    plt.plot(x, c, 'g.' , label='Temperature at each x')

    for t in range(1, 100):
        for i in range(1, n):
            a[i] = -(2 + (dx ** 2) )/ Dif(x[i]) * dt

            r[i] = -dx ** 2 * c[i] / Dif(x[i]) * dt
            b[i] = 1
            d[i] = 1
            delta[i] = - d[i] / (a[i] + b[i] * delta[i - 1])
            lambd[i] = (r[i] - b[i] * lambd[i - 1]) / (a[i] + b[i] * delta[i - 1])
        c[-1]=lambd[n]
    ii = n-1
    while(ii>0):
        c[ii]=delta[ii] * c[ii+1]+lambd[ii]
        ii=ii-1

    plt.plot(x, c, 'ro', label='Temperature at each xx')



    plt.title('Heat equation')
    plt.show()
    plt.legend()






