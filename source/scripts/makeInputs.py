import sys
import os

import numpy as np
import yaml

path = sys.argv[1]
outDir = sys.argv[2]
nParams = int(sys.argv[3])
nHidden = int(sys.argv[4])
miniMeth = sys.argv[5]
miniIter = int(sys.argv[6])
omegaList = map(float, sys.argv[7].strip('[]').split(',')) 
stepList = map(float, sys.argv[8].strip('[]').split(',')) 

for fname in os.listdir(path):
    with open(path + "/" + fname, 'r') as inputFile:
        tmpDict = yaml.load(inputFile, Loader=yaml.CLoader)
        omgIdx = omegaList.index(tmpDict["omega"])
        tmpDict["stepmc"] = stepList[omgIdx]
        tmpDict["numparameters"] = nParams
        tmpDict["numhiddenbias"] = nHidden
        tmpDict["maxitermc"] = 1048576
        tmpDict["importance"] = True
        tmpDict["jastrow"] = True
        tmpDict["progress"] = True
        tmpDict["minimization"] = [miniMeth, miniIter, 0.0001, 350000]
        tmpDict["resampling"] = ["autoblocking", 1]
        outFname = "w" + str(omegaList[omgIdx]) + "_D" + str(tmpDict["dim"]) +\
                "_N" + str(tmpDict["numparticles"]) + "_HFHO"
        tmpDict["output"] = outDir + "/" + outFname + ".txt"
        with open(outFname + "run.yaml", 'w') as outputFile:
            yaml.dump(tmpDict, outputFile)
