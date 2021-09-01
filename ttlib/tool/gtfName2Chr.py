import sys
sys.path.append('/home/xutun/software/ttlib')
from asEvents import gtfName2Chr
myArgs = sys.argv
gtfName2Chr(myArgs[1],myArgs[2],myArgs[3],int(myArgs[4]))