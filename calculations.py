from julia.api import Julia
from functions import *
from julia import Main

jpath = "C:/Users/adwait/AppData/Local/Programs/Julia-1.7.2/bin/julia.exe"
jl = Julia(runtime=jpath, compiled_modules=False) 

Main.include("utils.jl")
Main.include("wolff.jl")
'''
mag = jl.eval("wolff!(temp = 4.0, iters = 150)")

autocorr = auto_corr(np.array(mag), 100)

plt.plot(np.arange(len(autocorr)-1), autocorr[1:])
plt.grid()
plt.xlabel(r'$\Delta t$')
plt.ylabel(r'$\langle m(t) m(t+ \Delta t) \rangle - \langle m(t)^2 \rangle$')
plt.title('Auto-correlation function (T = 2) ')
plt.xticks(rotation= 40)
plt.show()
'''
temps, mags = jl.eval("diagram(wolff!)")

print(temps)
print(mags)