import matplotlib.pylab as plt
import pyalps,argparse
import os
import plottery

parser = argparse.ArgumentParser(description='Evaluate Renyi entropies for range of L', epilog='(C) Johannes Helmes 2015')

parser.add_argument('--infile','-i', help='Prefix of result files',required=True)
#parser.add_argument('--foreach','-f',default='h',help='Parameter name, (default h)')
parser.add_argument('--X','-x',default='h',help='Parameter name, (default h)')
parser.add_argument('--quantity','-q',default='h',help='Parameter name, (default IncNo)')
parser.add_argument('--plot','-p',action='store_true')
parser.add_argument('--reciprocal','-r',action='store_true')
parser.add_argument('--specheat','-s',action='store_true')
parser.add_argument('--verbose','-v',action='store_true')
parser.add_argument('--partitioned','-a',help="set this flag, if partitioned run",action='store_true')
args=parser.parse_args()

data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix=args.infile),[args.quantity,args.quantity+'2'])

if args.verbose:
    print data
plottery.empty_kiln()

dataG = pyalps.collectXY(data, x=args.X, y=args.quantity)
if args.specheat:
    y2=str(str(args.quantity)+'2')
    dataG2 = pyalps.collectXY(data, x=args.X, y=y2)
    b,e,yerr=plottery.convert_alps_dataset(dataG[0])
    b2,e2,yerr=plottery.convert_alps_dataset(dataG2[0])

    for beta,en,en2 in zip(b,e,e2):
        sp=(en2-(en**2))*beta**2
        print beta, sp

X=[]
Y=[]
Yerr=[]
for thing in dataG:
    x,y,yerr=plottery.convert_alps_dataset(thing)
    X.append(x)
    Y.append(y)
    Yerr.append(yerr)
    if not args.specheat:
        plottery.print_data(x,y,yerr)


if args.plot:
#plt.ylim([0, 2])
    plt.ylabel(args.quantity)
    plt.xlabel(args.X)
    plottery.alps_errorbar(dataG, fmt='o-')
    plottery.fire(type='sequential')
    plt.show()
