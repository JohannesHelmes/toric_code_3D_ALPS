import pyalps,argparse
import uncertainties as unc
import Evaluate_helper as eh
import numpy as np
import os


parser = argparse.ArgumentParser(description='Evaluate simple observables from ALPS simulations', epilog='(C) Johannes Helmes 2016')

parser.add_argument('--infile','-i', help='Prefix of result files',required=True)
#parser.add_argument('--foreach','-f',default='h',help='Parameter name, (default h)')
parser.add_argument('--X','-x',default='beta',help='Parameter name, (default beta)')
parser.add_argument('--quantity','-q',default='FullEnergy',help='Parameter name, (default FullEnergy)')
parser.add_argument('--N','-n',type=int,default=1000,help='Number of resamplings for bootstrap')
parser.add_argument('--reciprocal','-r',action='store_true')
parser.add_argument('--specheat','-s',action='store_true')
parser.add_argument('--magnsusc','-m',action='store_true')
parser.add_argument('--verbose','-v',action='store_true')
args=parser.parse_args()

def energy_flucs(e, e2, beta):
    return (e2 - e**2 ) * beta**2 

def magn_flucs(e, e2, beta):
    return (e2 - e**2 ) * beta

path=os.path.dirname(args.infile)
prefix=os.path.basename(args.infile)

N=args.N

dataG = eh.read_results_from_file(path, prefix, args.X, args.quantity)
if args.verbose:
    print dataG

if args.specheat:
    dataG2 = eh.read_results_from_file(path, prefix, args.X, args.quantity+'2')
    b, e, eerr = eh.convert_alps_dataset(dataG[0])
    b, e2, e2err= eh.convert_alps_dataset(dataG2[0])
    
    L = dataG[0].props['L']

    for beta,en,enerr,en2,en2err in zip(b,e,eerr,e2,e2err):
        cv, cverr = eh.bootstrap(energy_flucs, [en, en2, beta] , [enerr, en2err, 0.0], N=N)
        print beta, cv, cverr

elif args.magnsusc:
    dataG2 = eh.read_results_from_file(path, prefix, args.X, args.quantity+'2')
    b, e, eerr = eh.convert_alps_dataset(dataG[0])
    b, e2, e2err= eh.convert_alps_dataset(dataG2[0])
    
    L = dataG[0].props['L']

    for beta,en,enerr,en2,en2err in zip(b,e,eerr,e2,e2err):
        cv, cverr = eh.bootstrap(magn_flucs, [en, en2, beta] , [enerr, en2err, 0.0], N=N)
        print beta, cv, cverr

else:
    for thing in dataG:
        X,Y,Yerr=eh.convert_alps_dataset(thing)
        for (x,y,yerr) in zip(X,Y,Yerr):
            print x,y,yerr


