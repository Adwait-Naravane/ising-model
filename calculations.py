from functions import *

N = 100

T = np.arange(0.1, 5, 0.05)
config = latticegen(N)
sus = []
for t in T:
    susi = []
    for _ in range(10):
        config_fin = wolff(config, t, nsteps= 20)
        susi.append(np.var(config_fin)/t)
        config = latticegen(N)
    sus.append(np.mean(susi))

plt.scatter(T, sus, s = 10)
plt.plot(T, sus)
plt.grid()
plt.xlabel('Temperature')
plt.ylabel('Susceptibility')
plt.title('Wolff')
plt.show()