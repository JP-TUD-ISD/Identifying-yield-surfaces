import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import random
import os
import libmulti
import sys
import math
import optuna
from optuna.trial import TrialState
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from multiprocessing import Pool
from datetime import datetime
from pymoo.optimize import minimize
from pymoo.algorithms.soo.nonconvex.es import ES
from classforpymoo import VectorEvaluation
from pymoo.termination.fmin import MinimumFunctionValueTermination
from pymoo.termination import get_termination
from pymoo.core.termination import TerminateIfAny
from time import time
import copy as cp
import funfy

# feap interface
feapdir = ''  # path of FE software executable, if PythonFE is NOT usewd
# path to feap needs to be global path to work best
specifiers = ['s1', 's2', 's3']  # Parameters in the FEAP inputfile which need to get a sampled value
valtoreplace = 'u'  # identifier which is replaced by the sampled values
evaluationfile = 'fort.42169'  # file in which f_y is saved of the FEM software
histfile = 'fort.1478'  # file in which the history is stored
feapinputfile = 'i09.txt'  # name of the FEAP inputfile

# sampling in the eigenstress space
hydrorange = [-100, 100]  # interval of hydrostatic points
hydrosteps = 4  # how many points on the hydrostatic axis are tested
thetasteps = 100  # how many directions in the deviatoric plane are tested
namewalkdirections = 'wdir'  # file to store the walking directions
namehyddirectios = "hydro"  # Naming convention for files with hydrostatic information
hydroforprocs = True    # sampling strategy, if True: deviatoric directions are the same for each processor and hydrostatic points are distributed. If false, hydrostatic points are the same for all processors and deviatoric directions are distributed

# parameters for all optimizers
numericalzero = 0.1  # control if the found minimum is accepted as zero, optimization is aborted if that criterion is met
arange = [0, 200]  # the range in which the minima searched
anumbertrials = 12  # how many trials of Optuna are done per direction for optuna optimizer
alg = "nr"  # name of the algorithm, is overwritten by the if elif block at the end of the file
# for genetic algorithm, the above number is how many evaluations per generation are performed

# options for binary runs
dobinaryfiles = False  # use this with different FEAP input files, where the ifile stays and has an include to nameofbinary file
nameofbinary = "par.dat"

# pymoo genetic algorithm parameters
rule = 1.0 / 3.0    # factor of offsprings per generation, reducing the number of off-srpings in the first generation
n_max_gen=350        # maximum number of generations, after which the stop takes place

# pymoo general settings
generatehist = False  # Flag which can be set True or False, generates history for True, can be slower
verb = False  # set the verbosity of pymoo

# the arguments below can be used in a more sophisticated termination criterion
xtol=1e-20          # realtive tolerance in the x space
cvtol=1e-10         # cvtol
ftol=1e-8           # tolerance of the design space, realtive
period=2            # minimum number of periods, and amount of periods in ftol calc
n_max_evals=1000    # maximum number of function evaluations

# optuna parameters
initiala = 1.5  # start of the search
astep = 0.01  # step size in optuna
useastep = False

# NR Parameters
pert = 1.e-6  # distance between the two evaluations of fy, which are used to approximate the gradient via differences
countmax = 500  # maximum tries of Newtons's method before no identification is possible

# results
resultfilename = 'res.txt'  # name of the joint result file

# parallelization
numberproc = 2  # determines for how many processors the sampling is conducted and how many
# processes are started simultaneously
# is also used for HPC to determine the amount of cores the running scripts are created

# Python
virtualenvdir = 'venv'
usevirtualenv = True  # use python interpreter from virtualenvinronment, recommended!
suppressoutput = True  # FEAP output in terminal

# creating HPC script
taurusfilename = 'run-on-taurus.sh'
taurusruntime = '10:00:00'
partition = 'haswell64'
mailing = True
mailaddress = 'xx.yy@tu-dresden.de'
jobname = 'barnard-run'  # 'search-taurus' is alternative for running
# search taurus is compatible with Barnard cluster

# creating Barnard running script
barnardfilename = 'identi-on-barnard.sh'  # output file name
# runtime is taken from above
# mailing settings are used from above
# jobname is taken from above
barnardalgorithm = "es"  # alternatively optuna
# above defines which optimization strategy to use on the HPC

# Logfile settings changes
Lfilenamebase = 'Lfile'  # string for basis of log file
Lfileextension = '.log'  # extension of logfile
FlushLfile = True  # option if the log file is directly flushed after writing it, True better tracking of running jobs, but might be slower
Loggingonlyfailed = False

# HPC settings
usetempfolder = True  # only set this value true, if the HPC Cluster supports temp directories
useramfs = True  # this sets the directory which is used on the HPC Cluster, true means the local ramfs or tempfs directory is used
# always disable it for local runs
hpcname = "id-yield"  # assigning jobname on the new HPC

# plotting, deprecated, use pure-plot.py
planes = False

# anisotropic run
doanisotropic = False  # Flag to conduct search for minimum and maximum angle
rotationspecs = ['r1', 'r2', 'r3']  # values in the input, which are replaced by the angles
anglesteps = [30, 30, 30]  # steps in between minimum and maximum anlge
anglemin = [0.0, 0.0, 0.0]  # minimum value of sampled interval for angles, in degree
anglemax = [180.0, 180.0, 180.0]  # maximum value of sampled interval for rotation angles, degree

# settings for using FE Code with f2py
usepythonFE = True  # Flag to use the PythonFE code
pferotations = [37, 38, 39]  # the entries in material array, which denote rotation 
materialarray = [0.0,] * 500  # safety length
materialarray[0] = 6969.0  # matswitch
materialarray[6] = 5.5e5  # young's modulus
materialarray[7] = 0.4  # nu
materialarray[8] = 1.0  # idfy
materialarray[9] = 50.0  # Sigma0
materialarray[11] = 1.0  # idkap

# job to run, can be set on the terminal
do = 'plotting'
# 'debug-cut'
# 'optuna-search'
# 'sample-for-run'
# joinfiles
# plotting
# 'start-taurus'
# 'try-mp'
# 'testnumber'
# sample-for-run
# parallel-optuna
# ##########################################################################
# execution block
# for safety purposes
nodename = ""
corename = ""

# read in command line values from terminal
for i in range(1, len(sys.argv)):
    if i == 1:
        a = sys.argv[i]
        if not a.isnumeric():
            do = a
        else:
            print('No real running argument passed using: ' + do)
    elif i == 2:
        a = sys.argv[i]
        if a.isnumeric():
            numberproc = int(a)
    elif i == 3:
        a = sys.argv[i]
        if not a.isnumeric():
            hpcname = a
    elif i == 4:
        a = sys.argv[i]
        if not a.isnumeric():
            nodename = a
        else:
            nodename = str(a)
    elif i == 5:
        a = sys.argv[i]
        if not a.isnumeric():
            corename = a
        else:
            corename = str(a)

# constants they don't need to be adjusted except for debugging maybe
hydroaxis = np.ones(3)  # hydro static axis
startvec = np.zeros(3)  # initial vector for sampling in deviatoric plane, needs to be orthogonal to hydro static axis
startvec[0] = -1.0
startvec[2] = 1.0

# normalize vectors
hydronormed = hydroaxis / np.linalg.norm(hydroaxis)
startnormed = startvec / np.linalg.norm(startvec)

# initial values
hyd = 0.0
currhydropoint = [0.0, 0.0, 0.0]
walk = [0.0, 0.0, 0.0]

histrun = False
globaldir = " "

# time addition
tstart = time()
tlast = tstart
tnow = tstart

# logfile addition
logfilename = ""
flfile = ""


def timing():
    global tstart
    global tlast
    global tnow
    tnow = time()
    tinc = tnow - tlast
    tt = tnow - tstart
    tlast = tnow
    return tinc, tt


def logfiling(npro, behav, outforfile):
    global logfilename
    global Lfilenamebase
    global Lfileextension
    global flfile
    global FlushLfile
    global nodename
    global corename

    if behav == "init":
        logfilename = Lfilenamebase + str(npro)  # + Lfileextension
        if nodename:
            logfilename = logfilename + "--n-" + str(nodename)
        if corename:
            logfilename = logfilename + "--c-" + str(corename)
        logfilename = logfilename + Lfileextension
        flfile = open(file=logfilename, mode="w")
    elif behav == "log":
        tincrement, totaltime = timing()
        flfile.write(
            outforfile + "\t" + "time increment: " + str(tincrement) + "\t total time: " + str(totaltime) + "\n")
        if FlushLfile:
            flfile.flush()

    elif behav == "close":
        flfile.close()


# function to change the working directory to the /temp directory on the hpc, which limits
# usage of infinibandwith
def changedirhpc(direction):
    global globaldir
    global numberproc
    global histfile
    global feapinputfile
    global logfilename
    global useramfs
    global namehyddirectios
    global namewalkdirections
    global hydroforprocs

    if useramfs:
        di = "/dev/shm"
    else:
        di = "/tmp"
    newdir = di + "/" + str(hpcname) + "--" + str(
        numberproc) + "/"  # newdir = "/tmp/" + str(hpcname) + "--" + str(numberproc) + "/"
    testingdir = str(hpcname) + "--" + str(numberproc) + "/"
    if direction == 1:
        globaldir = os.getcwd()
        # testing if that directory exists
        if not os.path.exists(di):
            print("No /tmp directory found")
            exit()
        os.chdir(di)
        if not os.path.exists(str(testingdir)):
            os.mkdir(str(testingdir))
        os.chdir(globaldir)
        if hydroforprocs:
            os.system("cp " + str(feapinputfile) + " " + str(histfile) + " " + namewalkdirections + " " + "proc"
                      + str(numberproc) + namehyddirectios + " " + str(newdir))
        else:
            os.system("cp " + str(feapinputfile) + " " + str(histfile) + " " + namehyddirectios + " " + "proc"
                      + str(numberproc) + namewalkdirections + " " + str(newdir))
        os.chdir(newdir)
    if direction == 2 and globaldir == " ":
        print("the global directory was not defined when calling the change results back to global")
        return
    if direction == 2:
        os.system("cp *res" + " " + str(logfilename) + " " + str(globaldir))
        os.chdir(globaldir)


# rotate around an arbitrary axis, here used to rotate an start vector around hydro static axis
def rotaroundaxis(theta, rotaxis, vec):
    rotrot = np.zeros((3, 3))
    costheta = np.cos(theta)
    sintheta = np.sin(theta)
    ux = rotaxis[0]
    uy = rotaxis[1]
    uz = rotaxis[2]
    rotrot[0, 0] = costheta + ux * ux * (1 - costheta)
    rotrot[0, 1] = ux * uy * (1 - costheta) - uz * sintheta
    rotrot[0, 2] = ux * uz * (1 - costheta) + uy * sintheta
    rotrot[1, 0] = uy * ux * (1 - costheta) + uz * sintheta
    rotrot[1, 1] = costheta + uy * uy * (1 - costheta)
    rotrot[1, 2] = uy * uz * (1 - costheta) - ux * sintheta
    rotrot[2, 0] = uz * ux * (1 - costheta) - uy * sintheta
    rotrot[2, 1] = uz * uy * (1 - costheta) + ux * sintheta
    rotrot[2, 2] = costheta + uz * uz * (1 - costheta)
    return rotrot.dot(vec)


# set the correct output format for FEAP input files, since a limitation in characters exists
def format_e(n):
    asn = '%E' % n
    return asn.split('E')[0].rstrip('0').rstrip('.') + 'E' + asn.split('E')[1]


# function to replace the specifiers in the input file with sampled values
def changeinput(filenamefeap, spec, val, value):
    # file name is the name of the file to operate upon
    # spec is the list of specifiers, which indicate change
    # val the string, which is replaced
    # value contains the values, which are written as replacement, needs to be as long
    # as the spec list
    # alg from https://blog.finxter.com/how-to-search-and-replace-a-line-in-a-file-in-python/
    global dobinaryfiles
    global nameofbinary

    if not dobinaryfiles:
        infile = open(filenamefeap, "r")
        content = b''
        lines = infile.readlines()
        for line in lines:
            line = line.strip()
            identi = (line.split("=")[0]).strip()
            countl = 0
            for test in spec:
                if identi == test:
                    # print("found line for " + test)
                    var = cutinput(str(value[countl]))
                    line = line.replace(val, var)
                    break
                countl = countl + 1  # very inefficient, but works
            content = content + line.encode('utf-8') + b"\n"
        infile.close()
        infile = open(filenamefeap, "wb")
        infile.write(content)
        infile.close()
    else:
        infile = open(nameofbinary, "wb")
        infile.write(b"PARA\n")
        for i in range(len(spec)):
            infile.write(str(spec[i]).encode('utf-8') + b"=" + str(value[i]).encode("utf-8") + b"\n")
        infile.write(b"\n")
    return None

# new version of reducing the length of a variable, which will be written to a FEAP input file
def cutinput(inp):
    stringg = "%.8e" % float(inp)
    if 'e' in stringg:
        out = stringg.replace("e", "d")
    elif 'E' in stringg:
        out = stringg.replace("E", "d")
    else:
        out = stringg
    return out


# start FEAP from the terminal via os.system command
def startfeap(iname, feapname):
    global suppressoutput

    if suppressoutput:
        comm = feapname + " -i" + iname + " > /dev/null 2>&1"
    else:
        comm = feapname + " -i" + iname
    os.system(comm)

    return None


# function for evaluating a simulation
def evalfeap(evfile):
    calcs = np.loadtxt(fname=evfile)
    resres = calcs[0]
    kind = 'elastic'
    for res in calcs:
        if res > 0.0:
            kind = 'plastic'
        if res > resres:
            resres = res
    return kind, resres


# current version of a function that changes into the working directory, changes the input file, runs simulation and
# evaluates the simulation
def dofeap2(vall):
    global feapinputfile
    global www
    global specifiers
    global valtoreplace
    global feapdir
    global evaluationfile
    global histrun  # incorporated for follow up yield surfaces
    global dobinaryfiles

    if not dobinaryfiles:
        rempath = os.getcwd()
        os.system('cp ' + feapinputfile + ' "' + www + '"')
        if histrun:
            os.system('cp ' + histfile + ' "' + www + '"')
        os.chdir(www)
    changeinput(filenamefeap=feapinputfile, spec=specifiers, val=valtoreplace, value=vall)
    startfeap(iname=feapinputfile, feapname=feapdir)
    resi = evalfeap(evfile=evaluationfile)
    if not dobinaryfiles:
        os.system('rm ' + feapinputfile)
        if histrun:
            os.system('rm ' + histfile)
        os.chdir(rempath)
    return resi


# function which runs the optimization along one vector in the deviatoric plane based upon optuna
def objectiveschmective(trial):
    global arange
    global astep
    global feapinputfile
    global specifiers
    global valtoreplace
    global currhydropoint
    global walk
    global hyd
    global numericalzero
    global materialarray
    global usepythonFE

    if useastep:
        aaa = trial.suggest_float(name='aaa', low=arange[0], high=arange[1], step=astep)
    else:
        aaa = trial.suggest_float(name='aaa', low=arange[0], high=arange[1])

    # check if the trial state has been calculated already
    complete_trials = trial.study.get_trials(deepcopy=False, states=(TrialState.COMPLETE,))
    for t in complete_trials[::-1]:
        if trial.params == t.params:
            return t.value

    point = hyd * hydronormed + aaa * walk
    if not usepythonFE:
        behav, fy = dofeap2(vall=point)
    else:
        fy = funfy.getfy(materialarray, point[0], point[1], point[2], 0)
    if fy * fy < numericalzero:
        trial.study.stop()  # this let's the study close once a good approximation is found
        # the current value is still returned, so best result is included
    return fy * fy


def getfiles(pid):
    global filename1
    global www
    global usepythonFE
    global feapinputfile
    global histfile
    global histrun
    global hydroforprocs

    strnumberproc = str(pid)
    logfiling(pid, "init", "")
    if not usepythonFE:
        if not (os.path.isdir(strnumberproc)):
            os.mkdir(strnumberproc)
        os.chdir(strnumberproc)
        www = os.getcwd()
        if dobinaryfiles:
            os.chdir("..")
            os.system("cp " + str(feapinputfile) + " " + str(strnumberproc))
            if histrun:
                os.system("cp " + str(histfile) + " " + str(strnumberproc))
            os.chdir(strnumberproc)
            if hydroforprocs:
                filename1 = '../proc' + strnumberproc + namehyddirectios
                hydpoints = np.loadtxt(filename1)
                walkdirs = np.loadtxt("../" + namewalkdirections)
            else:
                filename1 = '../proc' + strnumberproc + namewalkdirections
                hydpoints = np.loadtxt("../" + namehyddirectios)
                walkdirs = np.loadtxt(filename1)
        else:
            if hydroforprocs:
                os.chdir("..")
                filename1 = 'proc' + strnumberproc + namehyddirectios
                hydpoints = np.loadtxt(filename1)
                walkdirs = np.loadtxt(namewalkdirections)
            else:
                os.chdir("..")
                filename1 = 'proc' + strnumberproc + namewalkdirections
                hydpoints = np.loadtxt(namehyddirectios)
                walkdirs = np.loadtxt(filename1)
    else:
        if hydroforprocs:
            filename1 = 'proc' + strnumberproc + namehyddirectios
            hydpoints = np.loadtxt(filename1)
            walkdirs = np.loadtxt(namewalkdirections)
        else:
            filename1 = 'proc' + strnumberproc + namewalkdirections
            hydpoints = np.loadtxt(namehyddirectios)
            walkdirs = np.loadtxt(filename1)
        www = os.getcwd()
    return hydpoints, walkdirs


def generalsearch(pid):
    global hyd
    global currhydropoint
    global walk
    global anumbertrials
    global numericalzero
    global dobinaryfiles
    global namewalkdirections
    global namehyddirectios
    global hydroforprocs
    global feapdir
    global feapinputfile
    global histfile
    global histrun
    global hydronormed
    global arange
    global evaluationfile
    global rule
    global xtol
    global cvtol
    global ftol
    global period
    global n_max_gen
    global n_max_evals
    global specifiers
    global valtoreplace
    global generatehist
    global verb
    global nameofbinary
    global materialarray
    global usepythonFE
    global countmax
    global initiala
    global pert
    global doanisotropic
    global alg

    if alg == "optuna":
        optuna.logging.set_verbosity(optuna.logging.WARNING)
    # strnumberproc = str(pid)
    resultsinfile = str(pid) + "res"
    # change for aniso, more res files
    resultfilenamemin = str(pid) + "resmin"
    resultfilenamemax = str(pid) + "resmax"

    hydpoints, walkdirs = getfiles(pid=pid)

    # catch the case of only one walkdirection
    if walkdirs.ndim == 1:
        walkdirs = np.reshape(walkdirs, (1, 3))
    # savepoints = []
    if hydpoints.ndim == 0:
        temp = np.zeros(1)
        temp[0] = hydpoints
        hydpoints = temp
    try:
        testing = hydpoints.shape[0]
    except IndexError:
        testing = 1
    totalsearches = testing * walkdirs.shape[0]
    count = 1

    if doanisotropic:
        savepointsmin, savepointsmax = anisoruns(pid=pid, walkdirs=walkdirs, hydpoints=hydpoints,
                                                 totalsearches=totalsearches, count=count)
        np.savetxt(fname=resultfilenamemin, X=savepointsmin)
        np.savetxt(fname=resultfilenamemax, X=savepointsmax)
    else:
        savepoints = doisoruns(pid=pid, walkdirs=walkdirs, hydpoints=hydpoints, totalsearches=totalsearches,
                               count=count)

        if savepoints:
            np.savetxt(fname=resultsinfile, X=savepoints)

    logfiling(pid, "close", "")
    return 0.0


def anglescalc():
    global anglesteps
    global anglemin
    global anglemax

    res = np.zeros((max(anglesteps), 3))

    for i in range(3):
        inc = (anglemax[i] - anglemin[i]) / (anglesteps[i] - 1)
        for j in range(anglesteps[i]):
            res[j, i] = j * inc

    return res


def anisoruns(pid, walkdirs, hydpoints, totalsearches, count):
    global hydroforprocs
    global hyd
    global currhydropoint
    global walk
    global anglesteps
    global usepythonFE
    global pferotations
    global Loggingonlyfailed

    savepointsmin = np.zeros((totalsearches, 3))
    savepointsmax = np.zeros((totalsearches, 3))

    cc = 0

    angles = anglescalc()

    if hydroforprocs:
        for hyd in hydpoints:
            currhydropoint = hyd * hydronormed
            fy = testhydro(ppoint=currhydropoint)  # test the hydropoint first
            if fy >= 0.0:  # good like this or numerical zero instead of 0.0
                print('hydropoint is plastic, saving only it')
                totalsearches = totalsearches - walkdirs.shape[0]
                for walk in walkdirs:
                    savepointsmin[cc, 0] = currhydropoint[0]
                    savepointsmin[cc, 1] = currhydropoint[1]
                    savepointsmin[cc, 2] = currhydropoint[2]
                    savepointsmax[cc, 0] = currhydropoint[0]
                    savepointsmax[cc, 1] = currhydropoint[1]
                    savepointsmax[cc, 2] = currhydropoint[2]
                    cc += 1
                if not Loggingonlyfailed:
                    lstring = "Hydrostatic point: " + str(hyd) + " plastic, saved hydrostatic point"
                    logfiling(pid, "log", lstring)
            else:
                for walk in walkdirs:
                    print('running study: ' + str(cc + 1) + 'from: ' + str(totalsearches))
                    # init min and max for each direction
                    aaamin = 1.0e300
                    aaamax = 0.0
                    for i in angles[:, 0]:
                        for j in angles[:, 1]:
                            for k in angles[:, 2]:
                                # set the angles
                                if usepythonFE:
                                    materialarray[pferotations[0]] = i
                                    materialarray[pferotations[1]] = j
                                    materialarray[pferotations[2]] = k
                                sp, aaa = runastudy(pid=pid, totalsearches=totalsearches, count=cc + 1)
                                if aaa > aaamax:
                                    aaamax = aaa
                                    savepointsmax[cc, 0] = sp[0]
                                    savepointsmax[cc, 1] = sp[1]
                                    savepointsmax[cc, 2] = sp[2]
                                elif aaa < aaamin:
                                    aaamin = aaa
                                    savepointsmin[cc, 0] = sp[0]
                                    savepointsmin[cc, 1] = sp[1]
                                    savepointsmin[cc, 2] = sp[2]
                    cc += 1
    else:
        for walk in walkdirs:
            for hyd in hydpoints:
                fy = testhydro(ppoint=currhydropoint)
                if fy >= 0.0:
                    print('hydropoint is plastic, saving only it')
                    totalsearches = totalsearches - 1
                    savepointsmin.append(currhydropoint)
                    savepointsmax.append(currhydropoint)
                    if not Loggingonlyfailed:
                        lstring = "Hydrostatic point: " + str(hyd) + " plastic, saved hydrostatic point"
                        logfiling(pid, "log", lstring)
                else:
                    print('running study: ' + str(count) + 'from: ' + str(totalsearches))
                    aaamin = 1.0e300
                    aaamax = 0.0
                    for i in angles[:, 0]:
                        for j in angles[:, 1]:
                            for k in angles[:, 2]:
                                # set the angles
                                if usepythonFE:
                                    materialarray[pferotations[0]] = i
                                    materialarray[pferotations[1]] = j
                                    materialarray[pferotations[2]] = k
                                sp, aaa = runastudy(pid=pid, totalsearches=totalsearches, count=cc + 1)
                                if aaa > aaamax:
                                    aaamax = aaa
                                    savepointsmax[cc, 0] = sp[0]
                                    savepointsmax[cc, 1] = sp[1]
                                    savepointsmax[cc, 2] = sp[2]
                                elif aaa < aaamin:
                                    aaamin = aaa
                                    savepointsmin[cc, 0] = sp[0]
                                    savepointsmin[cc, 1] = sp[1]
                                    savepointsmin[cc, 2] = sp[2]
                    cc += 1
    # print("gone here 2")
    return savepointsmin, savepointsmax


def doisoruns(pid, walkdirs, hydpoints, totalsearches, count):
    global hydroforprocs
    global hyd
    global currhydropoint
    global walk
    global Loggingonlyfailed

    savepoints = []

    cl = 1

    if hydroforprocs:
        for hyd in hydpoints:
            currhydropoint = hyd * hydronormed
            fy = testhydro(ppoint=currhydropoint)  # test the hydropoint first
            if fy >= 0.0:  # good like this or numerical zero instead of 0.0
                print('hydropoint is plastic, saving only it')
                totalsearches = totalsearches - walkdirs.shape[0]
                for walk in walkdirs:
                    savepoints.append(currhydropoint)
                if not Loggingonlyfailed:
                    lstring = "Hydrostatic point: " + str(hyd) + " plastic, saved hydrostatic point"
                    logfiling(pid, "log", lstring)
            else:
                for walk in walkdirs:
                    print('running study: ' + str(cl) + 'from: ' + str(totalsearches))
                    sp, aaa = runastudy(pid=pid, totalsearches=totalsearches, count=cl)
                    savepoints.append(sp)
                    cl += 1

    else:
        for walk in walkdirs:
            for hyd in hydpoints:
                fy = testhydro(ppoint=currhydropoint)
                if fy >= 0.0:
                    print('hydropoint is plastic, saving only it')
                    totalsearches = totalsearches - 1
                    savepoints.append(currhydropoint)
                    if not Loggingonlyfailed:
                        lstring = "Hydrostatic point: " + str(hyd) + " plastic, saved hydrostatic point"
                        logfiling(pid, "log", lstring)
                else:
                    print('running study: ' + str(cl) + 'from: ' + str(totalsearches))
                    sp, aaa = runastudy(pid=pid, totalsearches=totalsearches, count=cl)
                    savepoints.append(sp)
                    cl += 1
    return savepoints


def runastudy(pid, totalsearches, count):
    global filename1
    global www
    global hyd
    global currhydropoint
    global walk
    global anumbertrials
    global numericalzero
    global dobinaryfiles
    global namewalkdirections
    global namehyddirectios
    global hydroforprocs
    global feapdir
    global feapinputfile
    global histfile
    global histrun
    global hydronormed
    global arange
    global evaluationfile
    global rule
    global xtol
    global cvtol
    global ftol
    global period
    global n_max_gen
    global n_max_evals
    global specifiers
    global valtoreplace
    global generatehist
    global verb
    global nameofbinary
    global materialarray
    global usepythonFE
    global countmax
    global initiala
    global pert
    global alg
    global Loggingonlyfailed

    sp = []
    aout = 0.0

    if alg == "optuna":
        study = optuna.create_study(direction="minimize")
        study.optimize(objectiveschmective, n_trials=anumbertrials, gc_after_trial=True)
        founda = study.best_params.get('aaa')
        minimialfy = math.sqrt(study.best_value)
        if minimialfy <= numericalzero:
            tobe = hyd * hydronormed + founda * walk
            sp = tobe
            if not Loggingonlyfailed:
                lstring = 'study: ' + str(count) + 'from: ' + str(totalsearches) + " successful"
                logfiling(pid, "log", lstring)
            aout = founda
        else:
            lstring = 'study: ' + str(count) + 'from: ' + str(totalsearches) + " FAILED"
            logfiling(pid, "log", lstring)
            print("found no zero")
        count += 1
    elif alg == "ga":
        problem = VectorEvaluation(feapdir=feapdir, feapinputfile=feapinputfile,
                                   histfile=histfile, hyd=hyd,
                                   histrun=histrun, hydronormed=hydronormed, specifiers=specifiers,
                                   valtoreplace=valtoreplace,
                                   dobinaryfiles=dobinaryfiles, nameofbinary=nameofbinary, matarray=materialarray,
                                   usepythonFE=usepythonFE,
                                   www=www, walk=walk, xl=arange[0], xu=arange[1], evaluationfile=evaluationfile)

        if generatehist:
            algorithm = ES(n_offsprings=anumbertrials, rule=rule, save_history=True)
        else:
            algorithm = ES(n_offsprings=anumbertrials, rule=rule)
        # termination = DefaultSingleObjectiveTermination(xtol=xtol ,cvtol=cvtol ,ftol=ftol ,period=period ,n_max_gen=n_max_gen, n_max_evals=n_max_evals)
        terminator1 = get_termination("n_gen", n_max_gen)
        terminator2 = MinimumFunctionValueTermination(fmin=(numericalzero * numericalzero))
        arnold = TerminateIfAny(terminator1, terminator2)
        fy2 = minimize(problem=problem, algorithm=algorithm, termination=arnold, verbose=verb)
        founda = fy2.X
        minimialfy = math.sqrt(fy2.F)
        if generatehist:
            n_gen = len(fy2.history)
        if minimialfy <= numericalzero:
            tobe = hyd * hydronormed + cp.deepcopy(founda) * walk
            sp = tobe
            if not Loggingonlyfailed:
                if generatehist:
                   lstring = 'study: ' + str(count) + 'from: ' + str(totalsearches) + " successful, with : " + str(n_gen)  + " generations "
                else:
                   lstring = 'study: ' + str(count) + 'from: ' + str(totalsearches) + " successful"
                logfiling(pid, "log", lstring)
            aout = cp.deepcopy(founda)
        else:
            lstring = 'study: ' + str(count) + 'from: ' + str(totalsearches) + " FAILED"
            logfiling(pid, "log", lstring)
            print("found no zero")
        count += 1
        if generatehist:
            del fy2, arnold, algorithm, problem, terminator1, terminator2, minimialfy, founda, n_gen
        else:
            del fy2, arnold, algorithm, problem, terminator1, terminator2, minimialfy, founda
    elif alg == "nr":
        xn = initiala
        for i in range(countmax):
            aaa = xn
            point = hyd * hydronormed + aaa * walk
            if not usepythonFE:
                behav, fy = dofeap2(vall=point)
            else:
                fy = funfy.getfy(materialarray, point[0], point[1], point[2], 0)
            fofx = fy * fy
            if fofx < numericalzero:
                counts = i
                break
            aaap = xn + pert
            pointp = hyd * hydronormed + aaap * walk
            if not usepythonFE:
                behavp, fyp = dofeap2(vall=pointp)
            else:
                fyp = funfy.getfy(materialarray, pointp[0], pointp[1], pointp[2], 0)
            fofxp = fyp * fyp
            fprimeofx = (fofxp - fofx) / pert
            x1 = xn - fofx / fprimeofx
            xn = x1

        if fofx <= numericalzero:
            tobe = hyd * hydronormed + xn * walk
            sp = tobe
            if not Loggingonlyfailed:
                lstring = 'study: ' + str(count) + 'from: ' + str(totalsearches) + " successful, with: " + str(counts) + "tries "
                logfiling(pid, "log", lstring)
            aout = xn
        else:
            lstring = 'study: ' + str(count) + 'from: ' + str(totalsearches) + " FAILED"
            logfiling(pid, "log", lstring)
            print("found no zero")
        count += 1

    return sp, aout


def testhydro(ppoint):
    global usepythonFE
    global materialarray

    if not usepythonFE:
        behav, fy = dofeap2(vall=ppoint)
    else:
        fy = funfy.getfy(materialarray, ppoint[0], ppoint[1], ppoint[2], 0)
    return fy


# ##################################################################################
# checks if file etc. exists
if not usepythonFE:
    if not os.path.exists(feapinputfile):
        print("No valid input file")
        exit()
    if os.path.exists(histfile):
        print("run with initial history")
        histrun = True
# ##################################################################################

# large if conditions that determines the behaviour of the script
if do == 'test':  # plot deviatoric plane, used for: debugging of sampling
    hydronormed = hydroaxis / np.linalg.norm(hydroaxis)
    a = hydronormed.dot(startvec)
    print(a)
    if a != 0.0:
        print('Vectors are wrong')
        exit()
    else:
        print('Vectos are correct')

    # derived quantites
    thetainc = 2 * np.pi / thetasteps

    # testing the rotation matrix
    b = rotaroundaxis(theta=2 * np.pi, rotaxis=hydronormed, vec=startvec)
    print("initial vector")
    print(startvec)
    print('rotated around 2pi')
    print(b)
    b = rotaroundaxis(theta=np.pi, rotaxis=hydronormed, vec=startvec)
    print('rotated around pi')
    print(b)
    # testing the 3d plot
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    for i in range(thetasteps + 1):
        thet = i * thetainc
        testvec = rotaroundaxis(theta=thet, rotaxis=hydronormed, vec=startvec) + hydroaxis
        ax.scatter(testvec[0], testvec[1], testvec[2], marker='o')
    ax.scatter(hydroaxis[0], hydroaxis[1], hydroaxis[2], marker='o')
    plt.savefig("testplot.pdf")
    plt.show()

elif do == "testnumber":  # used for: tests of reducing length of input
    testnumber = random.random()
    print(testnumber)
    print(cutinput(testnumber))
    hihi = random.randint(-5, 5)
    testnumber = random.random() * pow(10, hihi)
    print(testnumber)
    print(cutinput(testnumber))

elif do == "run":  # test the rotation around an axis
    hydronormed = hydroaxis / np.linalg.norm(hydroaxis)
    startnormed = startvec / np.linalg.norm(startvec)
    thetainc = 2 * np.pi / thetasteps
    hydroinc = (hydrorange[1] - hydrorange[0]) / hydrosteps
    for i in range(hydrosteps + 1):
        currhydro = i * hydroinc + hydrorange[0]
        basevec = currhydro * hydronormed
        for ii in range(thetasteps + 1):
            currtheta = i * thetainc
            walkingdir = rotaroundaxis(theta=currtheta, rotaxis=hydronormed, vec=startnormed)
            # now get the foking function that searches along that dir
elif do == "try-mp":  # testing the shared memory parallelization interface
    # print(1)
    libmulti.CreateWorker(numberproc)

# Here the sampling of the input files takes place. This is done as a separate call of the script since in HPC
# environments, it makes more sense to do this in serial
elif do == "sample-for-run":  # Changed
    thetainc = 2 * np.pi / (thetasteps - 1)
    hydroinc = (hydrorange[1] - hydrorange[0]) / (hydrosteps - 1)
    savewalkdir = []
    for i in range(thetasteps):
        currtheta = i * thetainc  # theta starts at 0
        walkingdir = rotaroundaxis(theta=currtheta, rotaxis=hydronormed, vec=startnormed)
        savewalkdir.append(walkingdir)
    # since the same walking directions are tested for each point, they
    # only need to be saved in one file, reading is thread safe task later
    if hydroforprocs:
        np.savetxt(fname=namewalkdirections, X=savewalkdir)
    # sampling for the points on the hydrostatic axis
    savehydropoints = []
    for i in range(hydrosteps):
        currhydro = i * hydroinc + hydrorange[0]
        # currhydropoint = currhydro * hydronormed
        savehydropoints.append(currhydro)
    # save file per processor
    if hydroforprocs:
        # slice the array in sub arrays, sampled by multiprocessors
        elementsperprocessor = int(len(savehydropoints) / numberproc)
        for i in range(numberproc):
            filename = 'proc' + str(i + 1) + namehyddirectios
            tobesaved = savehydropoints[i * elementsperprocessor: (i + 1) * elementsperprocessor]
            np.savetxt(fname=filename, X=tobesaved)
    else:  # case of sampling the deviatoric vector per processork
        np.savetxt(fname=namehyddirectios, X=savehydropoints)
        elementsperprocessor = int(len(savewalkdir) / numberproc)
        for i in range(numberproc):
            filename = 'proc' + str(i + 1) + namewalkdirections
            tobesaved = savewalkdir[i * elementsperprocessor: (i + 1) * elementsperprocessor]
            np.savetxt(fname=filename, X=tobesaved)

# identification of yield surface for either application on HPC or serial on local machine with optuna optimization
elif do == "optuna-search":
    # debug = optunasearching(pid=numberproc)
    alg = "optuna"
    debug = generalsearch(pid=numberproc)

elif do == "ga-search":
    # debug = gasearch(pid=numberproc)
    alg = "ga"
    debug = generalsearch(pid=numberproc)

# run the identification of the yield function in shared memory parallelization
elif do == "parallel-optuna":
    alg = 'optuna'
    with Pool(numberproc) as pool:
        debug = pool.map(generalsearch, range(1, numberproc + 1))

# run the identification of the yield function in shared memory parallelization
elif do == "parallel-es":
    alg = "ga"
    with Pool(numberproc) as pool:
        debug = pool.map(generalsearch, range(1, numberproc + 1))

elif do == "nr-search":
    # debug = nrsearch(pid=numberproc)
    alg = "nr"
    debug = generalsearch(pid=numberproc)

elif do == "parallel-nr":
    alg = "nr"
    with Pool(numberproc) as pool:
        debug = pool.map(generalsearch, range(1, numberproc + 1))

# future feature: change on the hpc for the I/O heavy operations to a directory limiting network traffic
# NOT tested yet
# now barnard run supports two Optimization algorithms, optuna is random search, ES is evolutionary strategy from Pymoo
elif do == "barnard-run":
    if usetempfolder and not usepythonFE:
        changedirhpc(direction=1)
    if barnardalgorithm == "optuna":
        alg = "optuna"
        debug = generalsearch(pid=numberproc)
    elif barnardalgorithm == "es":
        alg = "ga"
        debug = generalsearch(pid=numberproc)
    elif barnardalgorithm == "nr":
        alg = "nr"
        debug = generalsearch(pid=numberproc)
    elif barnardalgorithm == "pes":
        alg = "ga"
        with Pool(numberproc) as pool:
            debug = pool.map(generalsearch, range(1, numberproc + 1))
    elif barnardalgorithm == "pnr":
        alg = "nr"
        with Pool(numberproc) as pool:
            debug = pool.map(generalsearch, range(1, numberproc + 1))
    else:
        print("No such optimization algorithm on barnard")
    if usetempfolder and not usepythonFE:
        changedirhpc(direction=2)

# join the results of each processor into one final result file
elif do == 'joinfiles':
    if not doanisotropic:
        if hydroforprocs:
            resfile = open(resultfilename, "w")
            for i in range(1, numberproc + 1):
                filetoopen = str(i) + "res"
                if not (os.path.exists(filetoopen)):
                    continue
                fff = open(filetoopen, "r")
                a = fff.readlines()
                for aa in a:
                    resfile.write(aa)
                fff.close()
        else:
            resfile = open(resultfilename, "w")
            elementsperprocessor = int(thetasteps / numberproc)
            for k in range(hydrosteps):
                for i in range(1, numberproc + 1):
                    # loop over all hydropoints
                    # now go to each processor and get all vectors belonging to that hydro static point
                    # which is from i * elementsperprocessor to (i+1) * elementsperprocessor
                    filetoopen = str(i) + "res"
                    if not (os.path.exists(filetoopen)):
                        continue
                    fff = open(filetoopen, "r")
                    a = fff.readlines()
                    for j in range(elementsperprocessor):
                        resfile.write(a[j * hydrosteps + k])
    else:
        if hydroforprocs:
            reee = resultfilename.split(".")[0]
            eeee = resultfilename.split(".")[-1]
            rmin = reee + "min." + eeee
            rmax = reee + "max." + eeee
            resfilemin = open(rmin, "w")
            for i in range(1, numberproc + 1):
                filetoopen = str(i) + "resmin"
                if not (os.path.exists(filetoopen)):
                    continue
                fff = open(filetoopen, "r")
                a = fff.readlines()
                for aa in a:
                    resfilemin.write(aa)
                fff.close()
            resfilemax = open(rmax, "w")
            for i in range(1, numberproc + 1):
                filetoopen = str(i) + "resmax"
                if not (os.path.exists(filetoopen)):
                    continue
                fff = open(filetoopen, "r")
                a = fff.readlines()
                for aa in a:
                    resfilemax.write(aa)
                fff.close()
        else:
            reee = resultfilename.split(".")[0]
            eeee = resultfilename.split(".")[-1]
            rmin = reee + "min." + eeee
            rmax = reee + "max." + eeee
            resfilemin = open(rmin, "w")
            elementsperprocessor = int(thetasteps / numberproc)
            for k in range(hydrosteps):
                for i in range(1, numberproc + 1):
                    # loop over all hydropoints
                    # now go to each processor and get all vectors belonging to that hydro static point
                    # which is from i * elementsperprocessor to (i+1) * elementsperprocessor
                    filetoopen = str(i) + "resmin"
                    if not (os.path.exists(filetoopen)):
                        continue
                    fff = open(filetoopen, "r")
                    a = fff.readlines()
                    for j in range(elementsperprocessor):
                        resfilemin.write(a[j * hydrosteps + k])
            resfilemax = open(rmax, "w")
            elementsperprocessor = int(thetasteps / numberproc)
            for k in range(hydrosteps):
                for i in range(1, numberproc + 1):
                    # loop over all hydropoints
                    # now go to each processor and get all vectors belonging to that hydro static point
                    # which is from i * elementsperprocessor to (i+1) * elementsperprocessor
                    filetoopen = str(i) + "resmax"
                    if not (os.path.exists(filetoopen)):
                        continue
                    fff = open(filetoopen, "r")
                    a = fff.readlines()
                    for j in range(elementsperprocessor):
                        resfilemax.write(a[j * hydrosteps + k])

# deprecated, create a shell script for old HPC
elif do == "start-taurus":
    scriptdir = os.getcwd() + "/mac-yield-func.py"
    tfile = open(taurusfilename, "w")
    tfile.write("#/bin/bash\n")
    tfile.write("#SBATCH --ntasks=" + str(numberproc) + "\n")
    tfile.write("#SBATCH --cpus-per-task=1\n")
    tfile.write("#SBATCH --mem-per-cpu=1G\n")
    tfile.write("#SBATCH --time=" + taurusruntime + "\n")
    tfile.write("#SBATCH --job-name=runpython\n")
    tfile.write("#SBATCH --partition=" + partition + "\n")
    if mailing:
        tfile.write("#SBATCH --mail-type=end\n")
        tfile.write("#SBATCH --mail-user=" + mailaddress + "\n")
    tfile.write("\n")
    tfile.write("scriptdir=" + scriptdir + "\n")
    tfile.write("\n")
    if usevirtualenv:
        tfile.write("module load Python/3.8.6\n")
        tfile.write("\n")
        tfile.write("source " + virtualenvdir + "/bin/activate \n")
        tfile.write("\n")
    else:
        tfile.write("module load modenv/hiera\n")
        tfile.write("module load GCCcore/10.2.0\n")
        tfile.write("module load Python/3.8.6\n")
        tfile.write("module load modenv/ml\n")
        tfile.write("module load SciPy-bundle/2022.05-foss-2022a\n")
        tfile.write("module load modenv/hiera\n")
        tfile.write("module load GCC/11.2.0\n")
        tfile.write("module load OpenMPI/4.1.1\n")
        tfile.write("module load matplotlib/3.4.3\n")
        tfile.write("\n")
    for i in range(1, numberproc + 1):
        tfile.write("srun --exclusive --ntasks=1 python $scriptdir " + jobname + " " + str(i) + " &\n")
        tfile.write("\n")
    tfile.write("\n")
    tfile.write('echo "Waiting for parallel jobs to complete"\n')
    tfile.write("wait \n")
    tfile.write('echo "All parallel jobs have been completed"\n')

# debugging the SLURM interface
elif do == "test-taurus":
    print(numberproc)

# debugging the SLURM interface on new HPC
elif do == "debug-barnard":
    print(numberproc)

# create a shell script to run the identification procedure on the new HPC
# the .sh file can be given to SLURM by sbatch <filename.sh>
elif do == 'create-input-barnard':
    scriptdir = os.getcwd() + "/mac-yield-func.py"
    tfile = open(barnardfilename, "w")
    tfile.write("#!/bin/bash\n")
    tfile.write("\n")
    tfile.write("#SBATCH --ntasks=" + "1" + "\n")  # only one due to array job
    tfile.write("#SBATCH --cpus-per-task=1\n")
    tfile.write("#SBATCH --mem-per-cpu=10G\n")
    tfile.write("#SBATCH --time=" + taurusruntime + "\n")
    tfile.write("#SBATCH --job-name=idyield-" + str(datetime.today().strftime('%Y-%m-%d')) + "\n")
    tfile.write("#SBATCH --array=1-" + str(numberproc) + "\n")  # run from 1 to n proc searches
    if mailing:
        tfile.write("#SBATCH --mail-type=end\n")
        tfile.write("#SBATCH --mail-user=" + mailaddress + "\n")
    tfile.write("\n")
    tfile.write("scriptdir=" + scriptdir + "\n")
    tfile.write("\n")
    if usevirtualenv:
        tfile.write("module purge\n")
        tfile.write("module load release/23.04\n")
        tfile.write("module load GCCcore/10.2.0\n")
        tfile.write("module load Python/3.8.6\n")
        tfile.write("\n")
        tfile.write("source " + virtualenvdir + "/bin/activate \n")
        tfile.write("\n")
    else:
        print("on barnard no configuration without virtual environments is programmed yet")
        exit()
    if usetempfolder:
        tfile.write("srun -n 1 python $scriptdir " + str(jobname) + " $SLURM_ARRAY_TASK_ID idyield-"
                    + str(datetime.today().strftime('%Y-%m-%d')) + " $SLURM_NODEID $SLURM_LOCALID \n")
    else:
        tfile.write(
            "srun -n 1 python $scriptdir " + str(jobname) + " $SLURM_ARRAY_TASK_ID $SLURM_NODEID $SLURM_LOCALID \n")
    tfile.write("\n")
    tfile.close()

# create a restart shell script
# can be used, if some processors did not complete the job due to wall time issues, then only the jobs that did not
# finish on the last run, are re-executed. Needs further debugging, since the array in #SBATCH can be too long for the line
elif do == "restart-barnard":
    # set all procs as active
    processors = []
    for i in range(1, numberproc + 1):
        processors.append(str(i))
    # filter out the processors, which produced an input file
    files = os.listdir()
    resfiles = []
    for i in files:
        if "res" in i and i != "res.txt":
            resfiles.append(i)

    for i in resfiles:
        number = int(i.split("res")[0])
        processors[number - 1] = "none"

    # create the final string, which should be written to the shell script
    string = ""
    finalstring = ""
    for count in range(len(processors)):
        i = processors[count]
        if i != "none":
            string = string + (i + ",")
        if count == numberproc - 1:
            finalstring = string[:-1]

    scriptdir = os.getcwd() + "/mac-yield-func.py"
    tfile = open(barnardfilename[:-3] + "rs.sh", "w")
    tfile.write("#!/bin/bash\n")
    tfile.write("\n")
    tfile.write("#SBATCH --ntasks=" + "1" + "\n")  # only one due to array job
    tfile.write("#SBATCH --cpus-per-task=1\n")
    tfile.write("#SBATCH --mem-per-cpu=10G\n")
    tfile.write("#SBATCH --time=" + taurusruntime + "\n")
    tfile.write("#SBATCH --job-name=idyield-" + str(datetime.today().strftime('%Y-%m-%d')) + "-rs\n")
    tfile.write("#SBATCH --array=" + finalstring + "\n")  # run from 1 to n proc searches
    if mailing:
        tfile.write("#SBATCH --mail-type=end\n")
        tfile.write("#SBATCH --mail-user=" + mailaddress + "\n")
    tfile.write("\n")
    tfile.write("scriptdir=" + scriptdir + "\n")
    tfile.write("\n")
    if usevirtualenv:
        tfile.write("module purge\n")
        tfile.write("module load release/23.04\n")
        tfile.write("module load GCCcore/10.2.0\n")
        tfile.write("module load Python/3.8.6\n")
        tfile.write("\n")
        tfile.write("source " + virtualenvdir + "/bin/activate \n")
        tfile.write("\n")
    else:
        print("on barnard no configuration without virtual environments is programmed yet")
        exit()
    if usetempfolder:
        tfile.write("srun -n 1 python $scriptdir " + str(jobname) + " $SLURM_ARRAY_TASK_ID idyield-"
                    + str(datetime.today().strftime('%Y-%m-%d')) + "\n")
    else:
        tfile.write("srun -n 1 python $scriptdir " + str(jobname) + " $SLURM_ARRAY_TASK_ID \n")
    tfile.write("\n")
    tfile.close()

# debug the reduction of length for the input file
elif do == "debug-cut":
    a = 50.0
    print(cutinput(a))

elif do == "check-out":
    log = ["none"] * numberproc
    for i in range(1, numberproc + 1):
        testfile = str(i) + "res"
        if os.path.isfile(testfile):
            f = open(testfile, "r")
            ln = f.readlines()
            lnlen = len(ln)
            f.close()
            if lnlen == (hydrosteps / numberproc * (thetasteps + 1)):
                log[i - 1] = "proc: " + str(i) + " passed \n"
            else:
                log[i - 1] = "proc: " + str(i) + " has only: " + str(lnlen) + " lines\n"
        else:
            log[i - 1] = "proc: " + str(i) + " no file\n"

    f = open("00_log.txt", "w")
    for i in log:
        f.write(i)
    f.close()

# deprecated plotting routine, use pure-plot.py
elif do == 'plotting':
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # ax.scatter(x, y, z, s=1, color='k')
    b = np.loadtxt(fname=resultfilename)
    x = b[:, 0]
    y = b[:, 1]
    z = b[:, 2]
    ax.scatter(x, y, z, s=1, color='k')
    ax.set_xlabel('sigma_1')
    ax.set_ylabel('sigma_2')
    ax.set_zlabel('sigma_3')
    # PolyPatches
    alpha = 0.2
    fc1 = "b"
    fc2 = "r"
    fc3 = "y"

    # x1 = [-0.6, -0.6, 0.15, 0.15]
    # y1 = [0.15, -0.6, -0.6, 0.15]
    # z1 = [0, 0, 0, 0]
    # z2 = [0.3, -0.3, -0.3, 0.3]

    x1 = [-100, -100, 100, 100]
    y1 = [100, -100, -100, 100]
    z1 = [0, 0, 0, 0]
    z2 = [100, -100, -100, 100]

    if planes:
        verts1 = [list(zip(x1, y1, z1))]
        verts2 = [list(zip(z1, x1, z2))]
        verts3 = [list(zip(x1, z1, z2))]

        pc1 = Poly3DCollection(verts1, edgecolors=fc1, linewidths=1)
        pc1.set_alpha(alpha)  # Order reversed
        pc1.set_facecolor(fc1)

        pc2 = Poly3DCollection(verts2, edgecolors=fc2, linewidths=1)
        pc2.set_alpha(alpha)  # Order reversed
        pc2.set_facecolor(fc2)

        pc3 = Poly3DCollection(verts3, edgecolors=fc3, linewidths=1)
        pc3.set_alpha(alpha)  # Order reversed
        pc3.set_facecolor(fc3)

        ax.add_collection3d(pc1)
        ax.add_collection3d(pc2)
        ax.add_collection3d(pc3)
    respdf = resultfilename.split(".")[0] + ".pdf"
    plt.savefig(respdf)
    plt.show()

# outline of the future features
# include a check, if the hydrostatic point itself already has f_y freater than zero
# increases calculation time
# -----------------------------------------------------------------------------------
# function that prints out the current state, to track progress, log file
