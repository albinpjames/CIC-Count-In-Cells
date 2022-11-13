import numpy as np
import math
import os


def testpp():
    from simulations import test_pp
    simul = test_pp()

    direc = simul.direc
    boxsize = simul.boxsize
    filename = simul.filename
    # Poisson Process
    import scipy.stats

    #Simulation window parameters
    xMin,xMax=0, boxsize
    yMin, yMax=0, boxsize
    zMin, zMax=0, boxsize
    xDelta, yDelta, zDelta =xMax-xMin, yMax-yMin, zMax-zMin  #rectangle dimensions
    areaTotal=xDelta*yDelta*zDelta
    
    #Point process parameters
    lambda0=1; #intensity (ie mean density) of the Poisson process
    
    #Simulate Poisson point process
    numbPoints = scipy.stats.poisson( lambda0*areaTotal ).rvs()#Poisson number of points
    xx = np.array(xDelta*scipy.stats.uniform.rvs(0,1,((numbPoints,1)))+xMin)#x coordinates of Poisson points
    yy = np.array(yDelta*scipy.stats.uniform.rvs(0,1,((numbPoints,1)))+yMin)#y coordinates of Poisson points
    zz = np.array(zDelta*scipy.stats.uniform.rvs(0,1,((numbPoints,1)))+zMin)

    mass = 10**(np.random.uniform(10, 15, size=(numbPoints)))
    print(mass)
    data = np.column_stack((mass,xx,yy,zz))

    directory = os.getcwd()
    pathset = os.path.join(directory,direc)

    # check the directory does not exist
    if not(os.path.exists(pathset)):
        os.makedirs(pathset)
    np.save(os.path.join(pathset,filename), data)

def test_rand():
    from simulations import test_rand
    simul = test_rand()

    direc = simul.direc
    boxsize = simul.boxsize
    filename = simul.filename
    
    #Simulate Poisson point process
    numbPoints = 1000000 #Number of points

    xx = np.random.uniform(0, boxsize, size=(numbPoints))
    yy = np.random.uniform(0, boxsize, size=(numbPoints))
    zz = np.random.uniform(0, boxsize, size=(numbPoints))

    mass = 10**(np.random.uniform(10, 15, size=(numbPoints)))
    print(mass)
    data = np.column_stack((mass,xx,yy,zz))

    directory = os.getcwd()
    pathset = os.path.join(directory,direc)

    # check the directory does not exist
    if not(os.path.exists(pathset)):
        os.makedirs(pathset)
    np.save(os.path.join(pathset,filename), data)

def abacuscosmos():
    from simulations import abacuscosmos
    import h5py

    simul = abacuscosmos()

    direc = simul.direc
    boxsize = simul.boxsize
    filename = simul.filename



    l_array=np.arange(0,24,1)
    mass, xx, yy, zz = [], [], [], []
    for l in l_array:
        fname = 'Data/emulator_1100box_planck_products/emulator_1100box_planck_spline_00_products/emulator_1100box_planck_spline_00_rockstar_halos/z1.000/halos/halos_0.'+str(l)+'.h5'
        f1=h5py.File(fname,'r+')
        x=f1['halos']
        pos=x['pos']
        mass.append(x['m'])
        xx.append(pos[:,0])
        yy.append(pos[:,1])
        zz.append(pos[:,2])

    data = np.column_stack((mass,xx,yy,zz))

    directory = os.getcwd()
    pathset = os.path.join(directory,direc)

    # check the directory does not exist
    if not(os.path.exists(pathset)):
        os.makedirs(pathset)
    np.save(os.path.join(pathset,filename), data)

def abacussummit():
    import abacussummit
    abacussummit.main()

if __name__ == "__main__":  
    #test_pp()
    test_rand()
    #abacussummit()
    #abacuscosmos()