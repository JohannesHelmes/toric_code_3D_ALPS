import matplotlib.pylab as plt
import pyalps,argparse
import os
import fnmatch
import numpy as np
from scipy.integrate import simps, romb
import Evaluate_helper as eh

parser = argparse.ArgumentParser(description='Perform thermodynamic integration to get O(1) contribution from Renyi entropies', epilog='(C) Johannes Helmes 2016')

parser.add_argument('--infile','-i', help='Prefix of result files',required=True)
parser.add_argument('--X','-x',default='h',help='Parameter name, (default h)')
parser.add_argument('--Y','-y',default='Magnetization',help='Parameter name, (default Magnetization)')
parser.add_argument('--Const','-c',default=np.log(2),type=float,help='Value of the constant for the integration, default=ln(2)')
parser.add_argument('--B','-b',default=1000,type=int,help='Number of Bootstrap samples, default=1000')
parser.add_argument('--N','-n',default=2,type=int,help='Renyi order, default=2')
parser.add_argument('--Dim','-d',default=2,type=int,help='Dimension of the system, default=2')
parser.add_argument('--System','-S',default='toriccode',choices=['toriccode','toriccode3D','ising','toriccodeupright'])
parser.add_argument('--verbose','-v',action='store_true')
parser.add_argument('--alternate','-a',action='store_true',help="If bipartitions 2 and 3 are NOT identical")
parser.add_argument('--simple','-s',type=int,help="With this flag, bipartition 0 is assumed to be empty, compute the entropy itself instead of gamma, provide the bipartition number > 0")
parser.add_argument('--Step','-e',default=2,type=int,help='Number of increment data points for thermodynamic integration, default=2')
parser.add_argument('--outfile','-o',help='Filename of the result file. If not given, print results to stdout.')

args=parser.parse_args()

def ReadResults(path, prefix, X, Y):
    """
    Input
    ------------------
    path:  path to result files
    prefix: prefix filenames (ending out.h5), if ending with a dot (.), only one simulation is assumed
                            else the separate simulation are identified and merged.
    X : string, name of parameter
    Y : string, name of observable

    returns
    -----------------
    pyalps X-Y-props list
    """

    if prefix[-1]=='.':
        dataset= pyalps.loadMeasurements(pyalps.getResultFiles(dirname=path,prefix=filename),args.Y)

        return pyalps.collectXY(dataset, x=X, y=Y, foreach=['IncNo'])
    else:
        all_prefixes=[]
        for f in os.listdir(path):
            if fnmatch.fnmatch(f, prefix+'*.out.h5'):
                before_first_dot=f.split('.')[0]+'.'
                if before_first_dot not in all_prefixes:
                    all_prefixes.append(before_first_dot)
        datasetsXY = []
        for pre in all_prefixes:
            tmp_dset= pyalps.loadMeasurements(pyalps.getResultFiles(dirname=path,prefix=pre),args.Y)
            datasetsXY.append( pyalps.collectXY(tmp_dset, x=X, y=Y, foreach=['IncNo']) )

        for i,dxy in enumerate(datasetsXY[0]):   #Take the 1st dataset as reference
            rIncNo = dxy.props['IncNo']
            for dslist in datasetsXY[1:]:
                for dxy2nd in dslist:
                    if rIncNo == dxy2nd.props['IncNo']:
                        datasetsXY[0][i] = pyalps.mergeDataSets([datasetsXY[0][i], dxy2nd])   

        return datasetsXY[0]

        


path=os.path.dirname(args.infile)
filename=os.path.basename(args.infile)
if args.verbose:
    print path, filename

if args.verbose:
    print pyalps.getResultFiles(prefix=args.infile)
    #print dataset

renyi_dataG = ReadResults(path, filename, args.X, args.Y)
n=args.N


if args.verbose:
    print renyi_dataG

if args.outfile!=None:
    outfile=open(args.outfile,'w')


factor={}
factor['ising']=2
factor['toriccode']=4
factor['toriccodeupright']=2
factor['toriccode3D']=3

N=args.B

def bootstrap(avg,err,N=1000):
    """ Performs a bootstrapping of an array of normal distributed data

    Parameters
    ----------
    avg : array_like
            Gaussian distributed average values
    err : array_like
            Errors on values
    N   : int, optional
            Number of bootstrapped values (default = 1000)

    Returns
    ---------
    samples : ndarray
            bootstrapped values, the result is an array of dimension avg.ndim + 1
    """

    #Needs numpy
    if len(avg) != len(err):
        raise ValueError("Arrays must be of the same dimension.")

    return np.array( [ np.random.normal(avg[j],err[j],size=N) for j in range(len(avg)) ] )


def thermodynamic_integration(x,y,yerr,N=1000):
    """ Performs the integration using Simpson's rule after bootstrapping the y values.
    """
    Q=np.zeros(N)

    y_bstrap=eh.resampling(y,yerr,N)
    for i in xrange(N):
        Q[i]=simps(y_bstrap[:,i],x,even='avg')

    return Q
    
    

def read_alps_data(alpsXY):
    """ Reads x, y, y_err from ALPS dataset and returns them
    """
    x=np.array(alpsXY.x)
    y_str=np.array(alpsXY.y)
    y=np.array([float(str(y_str[i]).split(" ")[0]) for i in range(len(y_str))])
    y_err=np.array([float(str(y_str[i]).split(" ")[2]) for i in range(len(y_str))])

    return x,y,y_err


NIncSteps=len(renyi_dataG)
S=np.zeros((NIncSteps,N))   #entropies for all bootstrapped samples
NValues=len(renyi_dataG[0].x)
L=renyi_dataG[0].props['L']

x=np.zeros((NIncSteps,NValues))
y=np.zeros((NIncSteps,NValues))
y_err=np.zeros((NIncSteps,NValues))

for schemes in renyi_dataG:
    scheme=int(schemes.props['IncNo'])
    x[scheme],y[scheme],y_err[scheme] = eh.convert_alps_dataset(schemes)


for k in range(3,NValues+1,args.Step):
    for schemes in renyi_dataG:
        scheme=int(schemes.props['IncNo'])
        S[scheme]=thermodynamic_integration(x[scheme,:k],y[scheme,:k],y_err[scheme,:k]+1e-15,N)

    if args.alternate:
        gamma=np.array([args.Const-(2*factor[args.System]* L ** args.Dim) * (S[0,l]-S[1,l]-S[2,l]+S[3,l]) for l in range(N)])
    elif args.simple:
        gamma=np.array([args.Const-(2*factor[args.System]* L ** args.Dim) * (S[args.simple,l]-S[0,l]) for l in range(N)])
    else:
        gamma=1./(n-1)*np.array([args.Const-(n*factor[args.System] * L ** args.Dim) * (S[0,l]-2*S[1,l]+S[2,l]) for l in range(N)])


    if args.outfile:
        outfile.write("{0:f} {1:f} {2:f} \n".format(x[0,k-1], np.mean(gamma), np.std(gamma)) )
    else:
        print x[0,k-1], np.mean(gamma), np.std(gamma)

            
if args.outfile!=None:
    outfile.close()
