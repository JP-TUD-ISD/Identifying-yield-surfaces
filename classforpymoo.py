import os
import numpy as np
import pymoo.optimize
from pymoo.core.problem import ElementwiseProblem
from pymoo.optimize import minimize
import shutil
from pymoo.algorithms.soo.nonconvex.es import ES
from pymoo.termination.default import DefaultSingleObjectiveTermination
import funfy

# class of pymoo package, which needs to run the base solution
class VectorEvaluation(ElementwiseProblem):
    def __init__(self, feapinputfile, feapdir, histfile, hyd, histrun, hydronormed, www, 
                 walk, specifiers, valtoreplace, dobinaryfiles, nameofbinary, matarray, usepythonFE,
                 xl:np.array, xu:np.array, n_var = 1, n_obj = 1,
                 evaluationfile: str = 'fort.42169', **kwargs):
        self.feapinputfile = feapinputfile
        self.feapdir = feapdir
        self.histfile = histfile
        self.hyd = hyd
        self.histrun = histrun
        self.hydronormed = hydronormed
        self.www = www
        self.walk = walk
        self.specifiers = specifiers
        self.valtoreplace = valtoreplace
        self.evaluationfile = evaluationfile
        self.dobinaryfiles = dobinaryfiles
        self.nameofbinary = nameofbinary
        self.suppressoutput = True
        self.matarray = matarray
        self.usepythonFE = usepythonFE
        super().__init__(n_var = n_var , n_obj = n_obj, xl = xl, xu = xu, **kwargs)
        
    def _evaluate(self, x, out, *args, **kwargs):
        point = self.hyd * self.hydronormed + x * self.walk
        if not self.usepythonFE:
            behav, fy = self.dofeap(vall=point)
        else:
            fy = funfy.getfy(self.matarray, point[0], point[1], point[2], 0)
        out['F'] = [fy**2]
        
    def dofeap(self, vall):
        # print(a)
        if not self.dobinaryfiles:
            rempath = os.getcwd()
            os.system('cp ' + self.feapinputfile + ' "' + self.www + '"')
            if self.histrun:
                os.system('cp ' + self.histfile + ' "' + self.www + '"')
            os.chdir(self.www)
        self.changeinput(vall=vall)
        self.startfeap()
        resi = self.evalfeap()
        if not self.dobinaryfiles:
            os.system('rm ' + self.feapinputfile)
            if self.histrun:
                os.system('rm ' + self.histfile)
            os.chdir(rempath)
        return resi
        
    def cutinput(self, inp):
        stringg = "%.8e" % float(inp)
        if 'e' in stringg:
            out = stringg.replace("e", "d")
        elif 'E' in stringg:
            out = stringg.replace("E", "d")
        else:
            out = stringg
        return out
    
    
    # start FEAP from the terminal via os.system command
    def startfeap(self):
        
        if self.suppressoutput:
            comm = self.feapdir + " -i" + self.feapinputfile + " > /dev/null 2>&1"
        else:
            comm = self.feapdir + " -i" + self.feapinputfile
        os.system(comm)
    
        return None
    
    
    # function for evaluating a simulation
    def evalfeap(self):
        calcs = np.loadtxt(fname=self.evaluationfile)
        resres = calcs[0]
        kind = 'elastic'
        for res in calcs:
            if res > 0.0:
                kind = 'plastic'
            if res > resres:
                resres = res
        return kind, resres
        
    def changeinput(self, vall):
        # file name is the name of the file to operate upon
        # spec is the list of specifiers, which indicate change
        # val the string, which is replaced
        # value contains the values, which are written as replacement, needs to be as long
        # as the spec list
        # alg from https://blog.finxter.com/how-to-search-and-replace-a-line-in-a-file-in-python/
        if not self.dobinaryfiles:
            infile = open(self.feapinputfile, "r")
            content = ''
            lines = infile.readlines()
            for line in lines:
                line = line.strip()
                identi = (line.split("=")[0]).strip()
                countl = 0
                for test in self.specifiers:
                    if identi == test:
                        # print("found line for " + test)
                        var = self.cutinput(str(vall[countl]))
                        line = line.replace(self.valtoreplace, var)
                        break
                    countl = countl + 1  # very inefficient, but works
                content = content + line + "\n"
            infile.close()
            infile = open(self.feapinputfile, "w")
            infile.write(content)
            infile.close()
        else:
            infile = open(self.nameofbinary, "wb")
            infile.write(b"PARA\n")
            for i in range(len(self.specifiers)):
                infile.write(str(self.specifiers[i]).encode('utf-8') + b"=" + str(vall[i]).encode("utf-8") + b"\n")
            infile.write(b"\n")
        return None

        
    
