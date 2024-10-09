#!/bin/bash

python3 -m numpy.f2py -c -m funfy -L/usr/lib/ -llapack assemblealgotangPplast.f calcfs.f calcPandMandel.f det3new.f dfinvdfsub.f dtrcbardfsub.f elment.f getfy.f Heavyside.f invert3new.f Piolamatlib.f plastpyieldfunctions.f plastwo.f rotateAB.f shpfunc.f stressandtangentisotropicneoP.f sweight.f 
python3 -m numpy.f2py assemblealgotangPplast.f calcfs.f calcPandMandel.f det3new.f dfinvdfsub.f dtrcbardfsub.f elment.f getfy.f Heavyside.f invert3new.f Piolamatlib.f plastpyieldfunctions.f plastwo.f rotateAB.f shpfunc.f stressandtangentisotropicneoP.f sweight.f -m funfy -h funfy.pyf
