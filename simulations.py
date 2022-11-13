from datetime import datetime

now = datetime.now()
timestamp = now.strftime("%d-%m-%Y_%H-%M-%S")
#timestamp = datetime.isoformat(now)

class test_pp:
    simname = "Test Data Poisson Process"
    '''Prameters from the simulation'''
    type, boxsize = 'test_pp', 20

    '''Filename to save'''
    filename_time = str(boxsize) + '_' + str(type)+ '_' + str(timestamp)
    filename = str(boxsize) + '_' + str(type)
    direc = "ProcessedData/TestPP"

    #cutsides = ['10','20','25']
    cutsides = ['10']

class test_rand:
    simname = "Test Data Poisson Process"
    '''Prameters from the simulation'''
    type, boxsize = 'test_rand', 100

    '''Filename to save'''
    filename_time = str(boxsize) + '_' + str(type)+ '_' + str(timestamp)
    filename = str(boxsize) + '_' + str(type)
    direc = "ProcessedData/TestRand"

    #cutsides = ['10','20','25']
    cutsides = ['10']

class abacussummit:
    simname = "Abaccus Summit"
    '''Prameters from the simulation'''
    #type, cosmo, intcont, boxsize, sf, ef = 'base', 'c000', 'ph000',2000, 0, 33
    type, cosmo, intcont, boxsize, sf, ef = 'small', 'c000', 'ph3000', 500, 0, 0

    '''Filename to save'''
    filename_time = str(boxsize) + "HMpc_" + cosmo + '_' + intcont + '_' + str(timestamp)
    filename = str(boxsize) + "HMpc_" + cosmo + '_' + intcont 
    direc = "ProcessedData/AbaccusSummit"

    #cutsides = ['10','20','25']
    cutsides = ['10']

    mass = 2.109081520453063*10**9

class abacuscosmos:
    simname = "Abaccus Cosmos"
    type, boxsize = 'AbcusCosmos', 720

    '''Filename to save'''
    filename_time = str(boxsize) + "HMpc_" + str(type)+ '_' + str(timestamp)
    filename = str(boxsize) + "HMpc_" + str(type) 
    direc = "ProcessedData/AbaccusCosmos"

    #cutsides = ['10','20','25']
    cutsides = ['10']

    mass = 2.109081520453063*10**9

