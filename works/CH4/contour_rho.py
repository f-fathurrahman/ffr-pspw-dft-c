import numpy as np
import matplotlib.pyplot as plt
import sys

# Set FFT grid size here (!HARDWIRED!)
NR1 = 26
NR2 = 22
NR3 = 22

data3d = open('STARTING_RHO.dat','r').readlines()

if(len(data3d) != NR1*NR2*NR3):
  print "Fourier grid size did not match with actual data"
  sys.exit(1)

rho1d = np.zeros(NR1*NR2*NR3)
for i in range(0,len(data3d)):
  rho1d[i] = float(data3d[i].split()[0])

# Reshape, using Fortran ordering
rho3d = np.reshape(rho1d, (NR1,NR2,NR3), order='F')

# Lattice parameters
a = 12.0
b = 10.0
c = 10.0
x = np.arange(0.0, a+0.1, a/(NR1-1))
y = np.arange(0.0, b+0.1, b/(NR2-1))
z = np.arange(0.0, c+0.1, c/(NR3-1))

# Plot in XY plane
X,Y = np.meshgrid(y,x) # the sequence is transposed?
plt.clf()
plt.figure(figsize=(6,6))
plt.contourf(X,Y,rho3d[:,:,0])
plt.savefig('contour1.png')


