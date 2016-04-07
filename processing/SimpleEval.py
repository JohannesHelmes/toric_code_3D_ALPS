import pyalps,argparse
import uncertainties as unc


parser = argparse.ArgumentParser(description='Evaluate simple observables from ALPS simulations', epilog='(C) Johannes Helmes 2016')

parser.add_argument('--infile','-i', help='Prefix of result files',required=True)
#parser.add_argument('--foreach','-f',default='h',help='Parameter name, (default h)')
parser.add_argument('--X','-x',default='beta',help='Parameter name, (default beta)')
parser.add_argument('--quantity','-q',default='FullEnergy',help='Parameter name, (default FullEnergy)')
parser.add_argument('--plot','-p',action='store_true')
parser.add_argument('--reciprocal','-r',action='store_true')
parser.add_argument('--specheat','-s',action='store_true')
parser.add_argument('--verbose','-v',action='store_true')
args=parser.parse_args()

def convert_alps_dataset(dset):
    X=dset.x
    Y = map(lambda x: float(str(x).split(" ")[0]), dset.y)
    Yerr = map(lambda x: float(str(x).split(" ")[2]), dset.y)
    return X,Y,Yerr

def convert_alps_dataset_unc(dset):
    X=dset.x
    Y = map(lambda x: unc.ufloat_fromstr(str(x)), dset.y)
    return X,Y

data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix=args.infile),[args.quantity,args.quantity+'2'])

if args.verbose:
    print data

dataG = pyalps.collectXY(data, x=args.X, y=args.quantity)
if args.specheat:
    y2=str(str(args.quantity)+'2')
    dataG2 = pyalps.collectXY(data, x=args.X, y=y2)
    #b,e,yerr=convert_alps_dataset(dataG[0])
    #b2,e2,yerr=convert_alps_dataset(dataG2[0])
    b, e = convert_alps_dataset_unc(dataG[0])
    b, e2, dump = convert_alps_dataset(dataG2[0])
    L = dataG[0].props['L']

    for beta,en,en2 in zip(b,e,e2):
        nquant= (en**2)
        print nquant.nominal_value, nquant.std_dev
        sp=(en2-(en**2))*beta**2
        print beta, sp.nominal_value, sp.std_dev

else:
    for thing in dataG:
        X,Y,Yerr=convert_alps_dataset(thing)
        for (x,y,yerr) in zip(X,Y,Yerr):
            print x,y,yerr


