#uses modified cclib 0.91
import sys
import logging 
from cclib.parser import Mopac
myfile=Mopac(sys.argv[1])
myfile.logger.setLevel(logging.ERROR) #cf. http://cclib.sourceforge.net/wiki/index.php/Using_cclib#Additional_information
data=myfile.parse()
print data.natom #print line 1
#print data.atomnos #print line 2; this has a problem with wrapping; Richard West suggested a couple of alternatives as shown below
print [i for i in data.atomnos] #this will put a comma between each element; alternative is:
#for i in data.atomnos:
#    print i,
print data.atomcoords[-1] #print the final optimized atom coordinates
print data.scfenergies[-1] #print the final optimized PM3 energy (originally was in eV, but I modified the modified cclib code so that it is in Hartree)
print data.molmass #print the molecular mass (in amu)
if (data.natom > 1):#for species with more than one atom, display vibrational frequencies and rotational constants
    #print data.vibfreqs #print the vibrational frequencies; see atomnos discussion above
    print [i for i in data.vibfreqs]  
    print data.rotcons[-1] #print the final rotational constants