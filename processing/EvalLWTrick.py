import matplotlib.pylab as plt
import pyalps,argparse
import numpy as np
import uncertainties
import uncertainties.umath
from collections import defaultdict

def renyi_sse_add(swap, for_each='L', **kwargs):
    if 'inc_range' not in kwargs:
        kwargs['inc_range'] = []

    values = defaultdict(dict)   

    for dset in swap:        
        if 'special' in kwargs:
            if kwargs['special'] == 'left_half':
                kwargs['inc_range'] = range(0, len(dset.x)/2);
            elif kwargs['special'] == 'right_half':
                kwargs['inc_range'] = range(len(dset.x)/2, len(dset.x));

        inc_range = kwargs['inc_range']

        if dset.props[for_each] in values.keys():
            pass
        else:
            values[dset.props[for_each]] = uncertainties.ufloat(0, 0)

        for x, y in zip(dset.x, dset.y):            
            if int(x) not in inc_range:
                continue

            u = uncertainties.ufloat_fromstr(str(y))                
            if u.nominal_value < 1e-12:
                print >> sys.stderr, "Plottery: Warning - entry", j, "is zero"
                continue
            values[dset.props[for_each]] += -uncertainties.umath.log(u/(uncertainties.ufloat_fromstr("1 +/- 0") -u))

    x_values = values.keys()
    y_values = [v.nominal_value for k, v in values.items()]
    y_errors = [v.std_dev for k, v in values.items()]

    return x_values, y_values, y_errors


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

renyi_dataG = pyalps.collectXY(data, x='IncNo', y='EG', foreach=[args.foreach])

if args.verbose:
    print renyi_dataG

if (args.steps!=None):
    IncNosIItoU=range(args.steps[0])
    IncNosUtoO=range(args.steps[0],args.steps[1])
    #IncNos = [%.1f % i for i in range(args.IncNoRange[0], args.IncNoRange[1])]
    print IncNosIItoU, IncNosUtoO
else:
    totIncNos=len(renyi_dataG[0].x)
    IncNosIItoU=range(totIncNos/2)
    IncNosUtoO=range(totIncNos/2,totIncNos)

X,Y1,Yerr1=renyi_sse_add(renyi_dataG, for_each=str(args.foreach),inc_name='IncNo',inc_range=IncNosIItoU)
X,Y2,Yerr2=renyi_sse_add(renyi_dataG, for_each=str(args.foreach),inc_name='IncNo',inc_range=IncNosUtoO)
Y=np.array(Y1)-np.array(Y2)
Yerr=np.array(Yerr1)+np.array(Yerr2)
for x,y,e in zip(X,Y,Yerr):
    print x, y, e

