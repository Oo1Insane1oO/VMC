import sys

def sn(n):
    """ return a string of n spaces """
    return "".join([" " for i in range(n)])
# end function sn

def tn(n):
    """ return a string of n tabs-spaces """
    return "".join(["    " for i in range(n)])
# end function tn

def allTrue(inputList):
    """ check if all elements in a list is True """
    t = False
    for l in inputList:
        if l:
            t = True
        else:
            t = False
            break
        # end ifelse
    # end forl
    return t
# end if

def recurrCheck(inputList, s, strings, i):
    """ check for string s recursively """
    if i >= len(inputList):
        return inputList, i
    # end ifeif

    if s in inputList[i-1]:
        """ join if right-adjacent """
        inputList[i-1] = " ".join(inputList[i-1:i+1])
        del inputList[i]
        if s in inputList[i-1].split()[-1]:
#         if allTrue([False if inputList[i-1].split()[-1] != j else True for j in
#             strings]):
            inputList, i = recurrCheck(inputList, s, strings, i)
        else:
            return inputList, i
    # end if
    return inputList , i
# end function recurrCheck

def conditionalConjoin(inputList, strings):
    """ join elements if elements are right after a string given in strings """
    i = 1;
    while (True):
        for string in strings:
            inputList, i = recurrCheck(inputList, string, strings, i)
        # end forS
        i += 1
        if i > len(inputList)-1:
            """ break condition """
            break
        # end if
    # end while True
    return inputList
# and conditionalConjoin

def replacer(inputString):
    """ remove redundant spaces """
    return inputString . replace("( ", "(", 10).replace(" )", ")").replace("("
            , "(", 10).replace(",  ", ", ")
# end function replacer

def readWrite(inputName, filename, className, numParameters):
    """ open scripts/templateWavefunction.string for reading and set filename.string """
    with open(inputName, "r") as readFile:
        """ read template file and set data """
        with open(filename, 'w') as openFile:
            """ open file for writing """
            for line in readFile:
                openFile.write(line.replace("claCL",
                    className).replace("clheader",
                        className.lower()).replace("NUMPART",
                            numParameters).replace("LCDEF",
                                className.upper()).replace("clabasis",
                                    className+"Basis"));
            # end for line
        # end with open openFile
    # end with open readFile

def makeHeaderSource(inputName, filename, className, numParameters):
    """ read template files, replace CL with input class name and write header
    and source """
    readWrite(inputName + ".h", filename + ".h", className, numParameters)
    readWrite(inputName + ".cpp", filename + ".cpp", className, numParameters)
# end function reader
 
try:
    path = sys.argv[1]
    className = sys.argv[2]
    numParameters = sys.argv[3]
#     jastrowFactor = sys.argv[4]
except IndexError:
#     print ("USAGE: python makeTemplateWavefunction.py 'path' 'classname'"
#             "'numParameters' 'jastrowFactor'")
    print ("USAGE: python makeTemplateWavefunction.py 'path' 'classname'"
            " 'numParameters'")
    sys.exit(1)
# end tryexcept

makeHeaderSource("scripts/templateWavefunction", path + ("/" if path[-1] != "/"
    else " ") + "wavefunctions/" + className.lower(), className, numParameters)
makeHeaderSource("scripts/templateBasis", path + ("/" if path[-1] != "/" else
    "") + "basis/" + className.lower() + "basis", className + "Basis",
    numParameters)
