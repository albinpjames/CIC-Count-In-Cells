import numpy as np
import matplotlib.pyplot as plt
import h5py

#halo mass is provided as number of particles
ParticleMass =     1e+10 #Msun h^-1
#for ploting mass fun
m_i=1.e10
m_f=1.e16  
bin_no=25
log_const=(np.log10(m_f)-np.log10(m_i))/bin_no
const=10**log_const
box=720.0001
n=5
w=n**3
k=box/n
#outbox=0
nullmass=0
hist=np.zeros(shape=(w,bin_no),dtype=float)
massfun=np.ndarray([])

halos = []
l_array=np.arange(0,32,1)
for l in l_array:
 fname = '/home/darkmatter/Desktop/emy/rockstar_halo/halos_0.'+str(l)+'.h5'
 f1=h5py.File(fname,'r+')
 x=f1['halos']
 pos=x['pos']
 mass=x['m']
 x1=pos[:,0]
 y1=pos[:,1]
 z1=pos[:,2]
 for ip,(d,e,f) in enumerate(zip(x1,y1,z1)):
  if m_i<mass[ip]<m_f :
   jx=int(d/k)
   jy=int(e/k)
   jz=int(f/k)
   b=int(jx+(n*jy)+(n*n*jz))
   p=int((np.log10(mass[ip])-np.log10(m_i))/log_const)
    #if p<bin_no and b<w:
   hist[b][p]=hist[b][p] + 1      
    #else :
     #outbox=outbox+1  
  else:
   nullmass=nullmass+1
  
print(hist)
print("nullmass",nullmass) 
b=np.empty((bin_no+1), dtype=float)
for i in range(0,bin_no+1):
 b[i]=(m_i)* (const)**i
#print(b)
#d is the bin width.
d=b[1:]-b[:-1]
#print('d=',d)
c=np.sqrt(b[1:]*b[:-1])
subvol=(box**3)/w
massfun=hist/(d*subvol)
print("massfun",massfun)
mean=np.mean(massfun,axis=0)
std_dev=np.std(massfun,axis=0)
print("mean",mean)
print("std_dev",std_dev)
A=(std_dev/mean)**2
#np.savetxt('mass',c)
np.savetxt('std_rock_halo5',A)
#f=np.sqrt(mean)
#g=f/(d*720**3.0)
#print(g)

#plt.errorbar(c,mean,yerr=g)
plt.errorbar(c,mean,yerr=std_dev,fmt='ro')
plt.yscale('log')
plt.xscale('log')
plt.show() 
