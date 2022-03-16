import numpy as np
import matplotlib.pyplot as plt
from numba import jit
import random



def latticegen(N):
  init_random = np.random.random((N,N))
  config = np.zeros((N,N))
  config[init_random >= 0.5] = 1
  config[init_random < 0.5] = -1
  return config

def wolff(config,T, nsteps=100):
  N = len(config)
  nbr = {(i,j):[((i+1)%N, j) , (i,(j+1)%N) , ((i-1)%N, j) , (i,(j-1)%N)] for i in range(N) for j in range(N)}
  p  = 1.0 - np.exp(-2.0 / T)
  for step in range(nsteps):
      k = np.random.randint(0,N)
      l = np.random.randint(0,N)
      Pocket, Cluster = [(k,l)], [(k,l)]
      while Pocket != []:
          j = random.choice(Pocket)
          for m in nbr[j]:
              if config[m] == config[j] and m not in Cluster \
                     and np.random.uniform(0.0, 1.0) < p:
                  Pocket.append(m)
                  Cluster.append(m)
          Pocket.remove(j)
      for j in Cluster:
          config[j] *= -1
  return config

@jit(nopython=True, nogil=True)
def mcmove(config, beta):
  N = len(config)
  i = np.random.randint(0, N)
  j = np.random.randint(0, N)
  s =  config[i, j]
  nb = config[(i+1)%N,j] + config[i,(j+1)%N] + config[(i-1)%N,j] + config[i,(j-1)%N]
  cost = 2*s*(nb)
  if cost < 0:
    s *= -1
  elif cost >= 0 and np.random.random() < np.exp(-cost*beta):
    s *= -1
  config[i, j] = s
  return config

def metropolis(config,beta, num_times = 10000):
  for _ in range(num_times):
    mcmove(config, beta)
  return config

def calcMag(config):
    #Magnetization of a given configuration
    mag = np.sum(config)
    return mag

def Energy_spin(config, i, j):
    S = config[i,j]
    nb = config[(i+1)%N, j] + config[i,(j+1)%N] + config[(i-1)%N, j] + config[i,(j-1)%N] 
    return -(nb)*S

def calcEnergy(config):
   #Energy of a given configuration (dimensionless)!
    energy = 0
    for i in range(len(config)):
        for j in range(len(config)):
            energy += Energy_spin(config, i, j)
    return energy/4

def auto_corr_fast(M, kappa):   
#   The autocorrelation has to be truncated at some point so there are enough
#   data points constructing each lag. Let kappa be the cutoff
    M = M - np.mean(M)
    N = len(M)
    fvi = np.fft.fft(M, n=2*N)
#   G is the autocorrelation curve
    G = np.real( np.fft.ifft( fvi * np.conjugate(fvi) )[:N] )
    G /= N - np.arange(N); G /= G[0]
    G = G[:kappa]
    return G