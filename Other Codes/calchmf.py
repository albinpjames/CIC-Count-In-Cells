import numpy as np
import math

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from scipy.stats import norm
from scipy.stats import poisson
import statistics

from calchmf import MassFunction
from abacusnbody.data.compaso_halo_catalog import CompaSOHaloCatalog
from astropy.cosmology import Planck18

from datetime import datetime

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
now = datetime.now()
timestamp = datetime.isoformat(now)

'''Redshifts to be computed'''
redshifts = ["3.000","0.200","0.100"]
#redshifts = ["8.000"]
#redshifts = ["0.200"]
#redshifts = ["8.000"]

'''Prameters from the simulation'''
type = 'base'
#type = 'small'
cosmo = 'c000'
intcont = 'ph000'
#intcont = 'ph3000'
boxsize = 2000
'''Filename to save'''
#filename = str(boxsize) + "HMpc " + cosmo + ' ' + intcont + ' ' + str(timestamp)
filename = str(boxsize) + ' ' + str(timestamp)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

'''Mass function from hmfcalc'''
def hmfcalc(redshift):
   print("HMF Calc") 
   global redshif
   redshif = float(redshift)
   global hmf
   hmf = MassFunction(cosmo_model=Planck18, z=redshift)
   return hmf  

#Read Data
def readdata(clm1, redshift):
   print("Reading the data")

   file = str('AbacusSummit Public Data Access/AbacusSummit_' + type + '_' + cosmo +'_' + intcont + '/halos/z' + redshift + '/halo_info/')
   cat = CompaSOHaloCatalog(file,cleaned=False)
   data = np.array(cat.halos[clm1])
   del cat
   return data 

def onebyone(clm1, redshift):
   print("Status: Reading the files")
   
   global halo_mass
   data = []
   for i in range(0,34): 
      print("Procesing file", i)
      if i < 10:
         file = str('AbacusSummit Public Data Access/AbacusSummit_' + type + '_' + cosmo +'_' + intcont + '/halos/z' + redshift + '/halo_info/halo_info_00' + str(i) + '.asdf')
      else:   
         file = str('AbacusSummit Public Data Access/AbacusSummit_' + type + '_' + cosmo +'_' + intcont + '/halos/z' + redshift +  '/halo_info/halo_info_0' + str(i) + '.asdf')
      cat = CompaSOHaloCatalog(file, cleaned=False)
      num_halos = len(cat.halos)
      num = np.array(cat.halos['N'])
      for i in range (num_halos):
         data.append(num[i])
      del cat
   
   return data

def convertmass(N):
   print("Converting mass")
   global num_halos, halo_mass
   halo_mass = []
   num_halos = len(N)
   for i in range (num_halos):
      halo_mass.append(N[i]*mass)

   total_num_halos = len(halo_mass)
   print("Number of halos =", total_num_halos)
   del N
   return halo_mass


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pdf = PdfPages('Plots/'+ filename + '.pdf')

fig = plt.figure(figsize=(20,30))
ax=[]
ax.append(plt.subplot2grid((3,2),loc=(0,0),rowspan=1,colspan=1))
ax.append(plt.subplot2grid((3,2),loc=(0,1),rowspan=1,colspan=1))
ax.append(plt.subplot2grid((3,2),loc=(1,0),rowspan=1,colspan=1))
ax.append(plt.subplot2grid((3,2),loc=(1,1),rowspan=1,colspan=1))
ax.append(plt.subplot2grid((3,2),loc=(2,0),rowspan=1,colspan=1))
ax.append(plt.subplot2grid((3,2),loc=(2,1),rowspan=1,colspan=1))

for j, redshift in enumerate(redshifts):
   print(redshift)
   
   mass = 2.109081520453063*10**9
   size = 2000
   #size = 500
   totalvolume = size**3

   totalbins = 100

   ax[0].axis('off')
   ax[0].text(0.05, 0.3, 'Simulation: AbacusSummit' +
                    '\nCosmology: Planck2018 LCDM' +
                    '\n\nBox Length: ' + str(size) + ' $\dfrac{MPc}{h}$' +
                    '\nParticle Mass:' + r'$2.1\times 10^9 \dfrac{M_\odot}{h}$' +
                    '\n\nNumber Of Mass Bins: ' + str(totalbins) 
                    ,size=22)


   #Extarct mass and position
   #N = readdata('N', redshift)
   N = onebyone('N', redshift)
   halo_mass = convertmass(N)
   del N

   binedge=np.logspace(np.log10(np.min(halo_mass))-0.5,np.log10(np.max(halo_mass))+0.5, totalbins+1)
   dm = np.log10(binedge[1])-np.log10(binedge[0])

   #Total HMF
   print("Calculating HMF")
   totalhmf = np.zeros(totalbins)
   for i in halo_mass:
      f = int((np.log10(i) - np.log10(binedge[0]))/dm)
      totalhmf[f]= totalhmf[f]+1

   totalbin_centers = []
   for i in range (len(binedge)-1):
      totalbin_centers.append(np.sqrt(binedge[i+1]*binedge[i]))

   global binsize 
   binsize = []
   for i in range (len(binedge)-1):
      binsize.append(binedge[i+1]-binedge[i])

   totalhmf_hght = []
   for i in range(len(totalhmf)):              
      totalhmf_hght.append(totalhmf[i]/((binsize[i]/totalbin_centers[i])*totalvolume))

   #nbin_centers, nhmf_hght = nbodyhmf(halo_mass,totalvolume)

   #Total HMF from hmfcalc
   hmf = hmfcalc(redshift)



   print('Plotting redshift:'+str(redshift))
   #Plotting

   #ax[j+1].scatter(nbin_centers, nhmf_hght, color = '#2ab0ff' , s=10)
   ax[j+1].scatter(totalbin_centers, totalhmf_hght, color = '#2ab0ff' , s=10, label='Computed HMF')
   ax[j+1].set_xscale('log')
   ax[j+1].set_yscale('log')
   xmin,xmax = ax[j+1].get_xlim()
   ymin,ymax = ax[j+1].get_ylim()
   ax[j+1].loglog(hmf.m,hmf.dndlnm, linewidth=0.5, color ='k', label='HMF Calc')
   ax[j+1].set_xlim(xmin,xmax)
   ax[j+1].set_ylim(ymin,ymax)
   ax[j+1].set_title("HMF for the redshift: " + str(redshift), fontsize=20)
   ax[j+1].set_xlabel("Mass $M$", fontsize=15)
   ax[j+1].set_ylabel("$\dfrac{dn}{Vdlnm}$", fontsize=15)
   ax[j+1].legend(loc="upper right")

   pdf.savefig(fig)
   #fig.savefig('Plots2/'+str(filename+'.png'))
   del halo_mass, totalbin_centers, totalhmf_hght
pdf.close()
