import sys
sys.path.append('/data/ttlib')
from searchNcbi import getAssemblyIdFromName
myArgs = sys.argv
print(getAssemblyIdFromName(myArgs[1]))