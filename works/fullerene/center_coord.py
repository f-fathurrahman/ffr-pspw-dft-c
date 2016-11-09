import numpy as np

geom = open('fullerene.xyz','r').readlines()

natom = int(geom[0].split()[0])
coord = np.zeros( (natom,3) )

atyp = []
ia=0
for il in range(2,natom+2):
  atyp.append(geom[il].split()[0])
  coord[ia,0] = float(geom[il].split()[1])
  coord[ia,1] = float(geom[il].split()[2])
  coord[ia,2] = float(geom[il].split()[3])
  ia = ia + 1

# Convert to bohr
coord = coord/0.529

# Dimension of box (in bohr)
a = 20.0
b = 20.0
c = 20.0
center = np.array([a/2.0, b/2.0, c/2.0])
cmx = (np.max(coord[:,0]) + np.min(coord[:,0]))/2.0
cmy = (np.max(coord[:,1]) + np.min(coord[:,1]))/2.0
cmz = (np.max(coord[:,2]) + np.min(coord[:,2]))/2.0

# Translate the coordinate
coord[:,0] = coord[:,0] + center[0]-cmx
coord[:,1] = coord[:,1] + center[1]-cmy 
coord[:,2] = coord[:,2] + center[2]-cmz

for ia in range(0,natom):
  print "%s %18.10f %18.10f %18.10f" % (atyp[ia], coord[ia,0], coord[ia,1], coord[ia,2])

