import sys

def degeneracy(i, dim):
    if (dim == 2):
        return 2*(i+1)
    else:
        return (i+1)*(i+2)
    # end ifelse
# end function degeneracy

def createMagic(numParticlesMax, dim):
    """ return list of magic number below incluse numParticlesMax """
    magic = [2]
    n = magic[0]
    while (n < numParticlesMax):
        magic.append(degeneracy(len(magic), dim))
        n += magic[-1]
    # end while
    for i in range(1,len(magic)):
        magic[i] += magic[i-1]
    # end fori
    return magic
# end function createMagic

def makeFiles(omegaList, numParticlesMax, dim, dName, dataDirName):
    particles = createMagic(numParticlesMax, dim)
    for w in omegaList:
        for p,i in enumerate(particles):
            fname = "w%.2f_D%i_N%i" % (w, dim, i)
            with open(dName + "/" + fname + "HOrun.yaml", "w") as openFile:
                openFile.write(('omega: %.2f\n'
                               'numparticles: %i\n'
                               'dim: %i\n'
                               'stepmc: 0.01\n'
                               'numparameters: 2\n'
                               'numhiddenbias: 0\n'
                               'maxitermc: 1048576\n'
                               'progress: true\n'
                               'importance: true\n' 
                               'jastrow: true\n'
                               'minimization: ["BFGS", 300, 0.00001, 215000]\n'
                               'resampling: ["autoblocking", 1]\n'
                               'output: '+'"' + dataDirName + '/' +
                               fname + ".yaml" + '"\n') % (w,i,dim))
            # end with open openFile
        # end fori
    # end forw
# end function makeFiles

if __name__ == "__main__":
    try:
        numParticlesMax = int(sys.argv[1])
        dim = int(sys.argv[2])
        dirName = sys.argv[3]
        dataDirName = sys.argv[4]
        omegaList = map(float, sys.argv[6:])
    except IndexError:
        print("USAGE: 'Nmax' 'dim' 'dir' 'datadir' 'list omega'")
        sys.exit(1)
    # end try-except

    makeFiles(omegaList, numParticlesMax, dim, dirName, dataDirName)
# end ifmain
