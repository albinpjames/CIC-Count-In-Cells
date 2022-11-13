import numpy as np
from abacusnbody.data.compaso_halo_catalog import CompaSOHaloCatalog
import os

from simulations import abacussummit
simul = abacussummit()
direc = simul.direc
filename = simul.filename

type, cosmo, intcont, boxsize, sf, ef = simul.type, simul.cosmo, simul.intcont, simul.boxsize, simul.sf, simul.ef

'''Redshifts to be computed'''
#redshifts = ["3.000","2.000","1.100","0.500","0.200"]
redshifts = ["3.000"]

mass = 2.109081520453063*10**9
#totalvolume = boxsize**3
totalbins = 50

def readdata(clm1, clm2, redshift):
   print("Reading the data")

   file = str('Data/AbacusSummit Public Data Access/AbacusSummit_' + type + '_' + cosmo +'_' + intcont + '/halos/z' + redshift + '/halo_info/')
   cat = CompaSOHaloCatalog(file,cleaned=False)
   data1 = np.array(cat.halos[clm1])
   data2 = np.array(cat.halos[clm2])
   del cat
   return data1, data2   
 

def onebyone(clm1, clm2, redshift):
   print("Status: Reading the files")
   
   global halo_mass
   data1 = []
   data2 = []
   for i in range(sf,ef+1): 
      print("Procesing file", i)
      if i < 10:
         file = str('AbacusSummit Public Data Access/AbacusSummit_' + type + '_' + cosmo +'_' + intcont + '/halos/z' + redshift + '/halo_info/halo_info_00' + str(i) + '.asdf')
      else:   
         file = str('AbacusSummit Public Data Access/AbacusSummit_' + type + '_' + cosmo +'_' + intcont + '/halos/z' + redshift +  '/halo_info/halo_info_0' + str(i) + '.asdf')
      cat = CompaSOHaloCatalog(file, cleaned=False)

      num_halos = len(cat.halos)
      num = np.array(cat.halos[clm1])
      for i in range (num_halos):
         data1.append(num[i])

      num_hpos = len(cat.halos)
      num = np.array(cat.halos[clm2])
      for i in range (num_hpos):
         data2.append(num[i])

      del cat
   
   return data1, data2

def convertmass(mass,N):
   print("Converting mass")
   global num_halos, halo_mass
   halo_mass = []
   num_halos = len(N)
   for i in range (num_halos):
      halo_mass.append(N[i]*mass)

   total_num_halos = len(halo_mass)
   print("Number of halos =", total_num_halos)
   del N
   halo_mass = np.array(halo_mass)
   halo_mass=halo_mass.reshape(halo_mass.shape[0],-1)
   return halo_mass

def main():
   for redshift in redshifts:
      print(redshift)
      #Extarct mass and position

      N, pos = onebyone('N','SO_central_particle', redshift)
      halo_mass = convertmass(mass,N)
      del N

      print("Extracting Position Info")
      halo_pos=[]
      for i in range(len(pos)):
         halo_pos.append(list(pos[i]))
      halo_pos=np.array(halo_pos)

      data = np.hstack((halo_mass,halo_pos))
      data[:,1:4] = data[:,1:4] + (boxsize/2)
      
      directory = os.getcwd()
      pathset = os.path.join(directory,"ProcessedData/AbaccusSummit",redshift)

      # check the directory does not exist
      if not(os.path.exists(pathset)):
         os.makedirs(pathset)
      print(filename)
      print(pathset)
      np.save(os.path.join(pathset,filename), data)

if __name__ == "__main__":  
   main()   