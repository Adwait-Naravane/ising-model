from julia.api import Julia
from functions import *
from julia import Main

jpath = "C:/Users/adwait/AppData/Local/Programs/Julia-1.7.2/bin/julia.exe"
jl = Julia(runtime=jpath, compiled_modules=False) 

Main.include("utils.jl")
Main.include("wolff.jl")

'''
mag = jl.eval("wolff!(temp = 2.0, iters = 270)")

autocorr = auto_corr_fast(np.array(mag), 700)

plt.plot(np.arange(len(autocorr)), autocorr)
plt.grid()
plt.xlabel(r'$\Delta t$')
plt.ylabel(r'$\langle m(t) m(t+ \Delta t) \rangle - \langle m(t)^2 \rangle$')
plt.title('Auto-correlation function (T = 2) ')
plt.xticks(rotation= 40)
plt.show()
'''
jl.eval("diagram(wolff!)")