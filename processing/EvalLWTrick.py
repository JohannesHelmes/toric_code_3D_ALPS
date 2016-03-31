import matplotlib.pylab as plt
import pyalps,argparse
import os
import plottery
import numpy as np

parser = argparse.ArgumentParser(description='Evaluate Renyi entropies for range of L', epilog='(C) Johannes Helmes 2014')

parser.add_argument('--infile','-i', help='Prefix of result files',required=True)
parser.add_argument('--foreach','-f',default='h',help='Parameter name, (default h)')
parser.add_argument('--steps','-s',nargs=2,type=int,help='Number of increment steps to complete U and O, (default=1/2 1)')
parser.add_argument('--plot','-p',action='store_true')
parser.add_argument('--verbose','-v',action='store_true')
args=parser.parse_args()

REntropy={}


data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix=args.infile),['EG'])


if args.verbose:
    print data
plottery.empty_kiln()

renyi_dataG = pyalps.collectXY(data, x='IncNo', y='EG', foreach=[args.foreach])

if args.verbose:
    print renyi_dataG

if (args.steps!=None):
    IncNosIItoU=range(args.steps[0])
    IncNosUtoO=range(args.steps[0],args.steps[1])
    #IncNos = [%.1f % i for i in range(args.IncNoRange[0], args.IncNoRange[1])]
else:
    totIncNos=len(renyi_dataG[0].x)
    IncNosIItoU=range(totIncNos/2)
    IncNosUtoO=range(totIncNos/2,totIncNos)

X,Y1,Yerr1=plottery.renyi_sse_add(renyi_dataG, for_each=str(args.foreach),inc_name='IncNo',inc_range=IncNosIItoU)
X,Y2,Yerr2=plottery.renyi_sse_add(renyi_dataG, for_each=str(args.foreach),inc_name='IncNo',inc_range=IncNosUtoO)
Y=np.array(Y1)-np.array(Y2)
Yerr=np.array(Yerr1)+np.array(Yerr2)
plottery.print_data(X,Y,Yerr)

if args.plot:
#plt.ylim([0, 2])
    plt.ylabel("$S_2(A)$")
    plt.xlabel(args.foreach)
    plottery.errorbar(X,Y,yerr=Yerr, fmt='o')
    plottery.fire(type='sequential')
    plt.show()
#plt.savefig("/home/helmes/Desktop/test_johanes.pdf")
