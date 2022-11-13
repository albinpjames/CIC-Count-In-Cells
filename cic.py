import numpy as np
import math

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from scipy.stats import norm
from scipy.stats import poisson
import statistics

from hmf import MassFunction
from abacusnbody.data.compaso_halo_catalog import CompaSOHaloCatalog
from astropy.cosmology import Planck18

from datetime import datetime

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
'''Redshifts to be computed'''
#redshifts = ["3.000","2.000","1.100","0.500","0.200"]
redshifts = ["3.000","5.000"]

'''Prameters from the simulation'''
#type, cosmo, intcont, boxsize, sf, ef = 'base', 'c000', 'ph000',2000, 0, 33
type, cosmo, intcont, boxsize, sf, ef = 'small', 'c000', 'ph3000', 500, 0, 0

now = datetime.now()
timestamp = now.strftime("%d/%m/%Y %H:%M:%S")
#timestamp = datetime.isoformat(now)

'''Filename to save'''
filename = str(boxsize) + "HMpc " + cosmo + ' ' + intcont + ' ' + str(timestamp)

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
def readdata(clm1, clm2, redshift):
   print("Reading the data")

   file = str('AbacusSummit Public Data Access/AbacusSummit_' + type + '_' + cosmo +'_' + intcont + '/halos/z' + redshift + '/halo_info/')
   cat = CompaSOHaloCatalog(file,cleaned=False)
   data1 = np.array(cat.halos[clm1])
   data2 = np.array(cat.halos[clm2])
   del cat
   return data1, data2   
 

def onebyone(clm1, clm2,  redshift):
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
   return halo_mass

#Dividing Boxes
def calcboxes (size, cutside, halo_pos):
   print("Divideing boxes: "+ str(cutside))

   nw_size = size/cutside
   x=(halo_pos + (size/2))/nw_size
   x = x.astype(int)
   boxes = []
   for i in range (len(x)):
      boxes.append(x[i][0] + cutside*x[i][1] + cutside**2*x[i][2])
   del halo_pos
   return boxes

def nbodyhmf(halo_mass,volume):
   print("Status: Hmf")
   print("Status: Calculating height and width of mass histogarm")
   hght, binedge= np.histogram(halo_mass, bins=np.logspace(np.log10(10**10),np.log10(10**15), 101,))

   print("Status: Calculating the mass function")
   bin_centers = []
   for i in range (len(binedge)-1):
      bin_centers.append(np.sqrt(binedge[i+1]*binedge[i]))
   
   global binsize 
   binsize = []
   for i in range (len(binedge)-1):
      binsize.append(binedge[i+1]-binedge[i])
 
   hmf_hght = []
   for i in range(len(hght)):              
      hmf_hght.append(hght[i]/((binsize[i]/bin_centers[i])*volume))
   return bin_centers, hmf_hght

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

for redshift in redshifts:
   print(redshift)
   #size, cutsides = 2000, ['40','80','100']
   size, cutsides = 500, ['10','20','25']

   mass = 2.109081520453063*10**9
   totalvolume = size**3
   totalbins = 50

   #Extarct mass and position

   N, pos = onebyone('N','SO_central_particle', redshift)
   halo_mass = convertmass(mass,N)
   del N

   print("Extracting Position Info")
   halo_pos=[]
   for i in range(len(pos)):
      halo_pos.append(list(pos[i]))
   halo_pos=np.array(halo_pos)

   #pdf = PdfPages('Plots/'+ str(redshift) +'z ' + filename + '.pdf')
   pdf = PdfPages(str(redshift) +'z ' + filename + '.pdf')

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

   for cutside in cutsides:
      cutside = int(cutside)
      totalboxes = cutside**3
      boxvolume = totalvolume/totalboxes

      boxes = calcboxes(size, cutside,halo_pos)

      #Box halo numbers
      data = np.zeros((totalbins,totalboxes))
      for i in range(len(boxes)):
         f = int((np.log10(halo_mass[i]) - np.log10(binedge[0]))/dm)
         data[f][boxes[i]] = data[f][boxes[i]]+1

      print('Plotting cutside:' + str(cutside))
      #Plotting
      fig = plt.figure(figsize=(20,30))
      ax0=plt.subplot2grid((3,2),loc=(0,0),rowspan=1,colspan=1)
      ax1=plt.subplot2grid((3,2),loc=(0,1),rowspan=1,colspan=1)
      ax2=plt.subplot2grid((3,2),loc=(1,0),rowspan=1,colspan=1)
      ax3=plt.subplot2grid((3,2),loc=(1,1),rowspan=1,colspan=1)
      ax4=plt.subplot2grid((3,2),loc=(2,0),rowspan=1,colspan=1)
      ax5=plt.subplot2grid((3,2),loc=(2,1),rowspan=1,colspan=1)


      #Plot HMF
      #ax1.scatter(nbin_centers, nhmf_hght, color = '#2ab0ff' , s=10)
      ax1.scatter(totalbin_centers, totalhmf_hght, color = '#2ab0ff' , s=10)
      ax1.set_xscale('log')
      ax1.set_yscale('log')
      ax1.axvline(x=10**11, color='#2ab0ff')
      ax1.axvline(x=5*10**11, color='#2ab0ff')
      ax1.axvline(x=10**12, color='#2ab0ff')
      ax1.axvline(x=5*10**12, color='#2ab0ff')
      xmin,xmax = ax1.get_xlim()
      ymin,ymax = ax1.get_ylim()
      ax1.loglog(hmf.m,hmf.dndlnm, linewidth=0.5, color ='k')
      ax1.set_xlim(xmin,xmax)
      ax1.set_ylim(ymin,ymax)
      ax1.set_title("HMF of the total box", fontsize=20)
      ax1.set_xlabel("Mass $M$", fontsize=15)
      ax1.set_ylabel("$\dfrac{dn}{Vdlnm}$", fontsize=15)

      bins = 'sturges'

      #Plot for mass 10**11
      mass = 10**11
      f=int((np.log10(mass) - np.log10(binedge[0]))/dm)
      ax2.hist(data[f][:], bins=bins, density= True, facecolor = '#2ab0ff', edgecolor='#169acf', linewidth=0.5)
      ax2.set_title(f"{mass:.1E}" +" (bin " + str(f)+")  Bins in this plot:"+str(bins), fontsize=20)
      ax2.set_xlabel("Number of halos in a box", fontsize=15)
      ax2.set_ylabel("Normalised Frequency", fontsize=15)

      xmin, xmax = ax2.get_xlim()
      xmin=0
      x_box_no = np.arange(xmin, xmax, 0.01)
      x_poi_no = np.arange(xmin, xmax)
      mean = statistics.mean(data[f][:])
      sd = statistics.stdev(data[f][:])
      ax2.plot(x_box_no, norm.pdf(x_box_no, mean, sd), 'k--', label='Gaussian')
      ax2.plot(x_poi_no, poisson.pmf(x_poi_no, mean), 'k', label='Poisson' )
      ax2.legend(loc="upper right")
      

      #Plot for mass 5*10**11
      mass = 5*10**11
      f=int((np.log10(mass) - np.log10(binedge[0]))/dm)
      ax3.hist(data[f][:], bins=bins, density= True, facecolor = '#2ab0ff', edgecolor='#169acf', linewidth=0.5)
      ax3.set_title(f"{mass:.1E}" +" (bin " + str(f)+")  Bins in this plot:"+str(bins), fontsize=20)
      ax3.set_xlabel("Number of halos in a box", fontsize=15)
      ax3.set_ylabel("Normalised Frequency", fontsize=15)

      xmin, xmax = ax3.get_xlim()
      xmin=0
      x_box_no = np.arange(xmin, xmax, 0.01)
      x_poi_no = np.arange(xmin, xmax)
      mean = statistics.mean(data[f][:])
      sd = statistics.stdev(data[f][:])
      ax3.plot(x_box_no, norm.pdf(x_box_no, mean, sd), 'k--', label='Gaussian')
      ax3.plot(x_poi_no, poisson.pmf(x_poi_no, mean), 'k', label='Poisson' )
      ax3.legend(loc="upper right")


      #Plot for mass 10**12
      mass =10**12
      f=int((np.log10(mass) - np.log10(binedge[0]))/dm)
      ax4.hist(data[f][:], bins=bins, density= True, facecolor = '#2ab0ff', edgecolor='#169acf', linewidth=0.5)
      ax4.set_title(f"{mass:.1E}" +" (bin " + str(f)+")  Bins in this plot:"+str(bins), fontsize=20)
      ax4.set_xlabel("Number of halos in a box", fontsize=15)
      ax4.set_ylabel("Normalised Frequency", fontsize=15)

      xmin, xmax = ax4.get_xlim()
      xmin=0
      x_box_no = np.arange(xmin, xmax, 0.01)
      x_poi_no = np.arange(xmin, xmax)
      mean = statistics.mean(data[f][:])
      sd = statistics.stdev(data[f][:])
      ax4.plot(x_box_no, norm.pdf(x_box_no, mean, sd), 'k--', label='Gaussian')
      ax4.plot(x_poi_no, poisson.pmf(x_poi_no, mean), 'k', label='Poisson' )
      ax4.legend(loc="upper right")

      #Plot for mass 5*10**12
      mass = 5*10**12
      f=int((np.log10(mass) - np.log10(binedge[0]))/dm)
      ax5.hist(data[f][:], bins=bins, density= True, facecolor = '#2ab0ff', edgecolor='#169acf', linewidth=0.5)
      ax5.set_title(f"{mass:.1E}" +" (bin " + str(f)+")  Bins in this plot:"+str(bins), fontsize=20)
      ax5.set_xlabel("Number of halos in a box", fontsize=15)
      ax5.set_ylabel("Normalised Frequency", fontsize=15)

      xmin, xmax = ax5.get_xlim()
      xmin=0
      x_box_no = np.arange(xmin, xmax, 0.01)
      x_poi_no = np.arange(xmin, xmax)
      mean = statistics.mean(data[f][:])
      sd = statistics.stdev(data[f][:])
      ax5.plot(x_box_no, norm.pdf(x_box_no, mean, sd), 'k--', label='Gaussian')
      ax5.plot(x_poi_no, poisson.pmf(x_poi_no, mean), 'k', label='Poisson' )
      ax5.legend(loc="upper right")

      
      ax0.axis('off')
      ax0.text(0.05, 0.3, 'Simulation: AbacusSummit' +
                           '\nCosmology: LCDM' +
                           '\nRedshift: ' + str(redshift) + 
                           '\nBox Length: ' + str(size) + ' MPc/h' +
                           '\n\nSub Box Length: ' + str(size/cutside) + ' MPc/h' +
                           '\nTotal Number Of Sub Boxes: ' + str(totalboxes) +
                           '\n\nNumber Of Mass Bins: ' + str(totalbins) +
                           '\n\nNumber of halos: ' + str(len(halo_mass))
                           ,size=22)

      pdf.savefig(fig)
   del halo_mass, halo_pos, data, totalbin_centers, totalhmf_hght
   pdf.close()
