import matplotlib.pylab as plt
import pyalps,argparse
import os
import plottery
import numpy as np
from scipy.integrate import simps, romb

parser = argparse.ArgumentParser(description='Perform thermodynamic integration to get O(1) contribution from Renyi entropies', epilog='(C) Johannes Helmes 2014')

parser.add_argument('--infile','-i', help='Prefix of result files',required=True)
parser.add_argument('--X','-x',default='h',help='Parameter name, (default h)')
parser.add_argument('--Y','-y',default='Magnetization',help='Parameter name, (default Magnetization)')
parser.add_argument('--Const','-c',default=np.log(2),type=float,help='Value of the constant for the integration, default=ln(2)')
parser.add_argument('--B','-b',default=1000,type=int,help='Number of Bootstrap samples, default=1000')
parser.add_argument('--N','-n',default=2,type=int,help='Renyi order, default=2')
parser.add_argument('--Dim','-d',default=2,type=int,help='Dimension of the system, default=2')
parser.add_argument('--System','-s',default='toriccode',help="['toriccode','toriccode3D','ising','toriccodeupright']")
parser.add_argument('--verbose','-v',action='store_true')
parser.add_argument('--Simple','-S',action='store_true')
parser.add_argument('--alternate','-a',action='store_true',help='Use this flag, if bipartition 2 and 3 are different, and 4 bipartitions are given')
parser.add_argument('--Step','-e',default=2,type=int,help='Number of increment data points for thermodynamic integration, default=2')
parser.add_argument('--outfile','-o',help='Filename of the result file. If not given, print results to stdout.')

args=parser.parse_args()

REntropy={}

path=os.path.dirname(args.infile)
filename=os.path.basename(args.infile)
print path, filename
dataset= pyalps.loadMeasurements(pyalps.getResultFiles(dirname=path,prefix=filename),args.Y)

if args.verbose:
    print pyalps.getResultFiles(prefix=args.infile)
    #print dataset

renyi_dataG = pyalps.collectXY(dataset, x=args.X, y=args.Y, foreach=['IncNo'])
n=args.N


if args.verbose:
    print renyi_dataG

factor={}
factor['ising']=2
factor['toriccode']=4
factor['toriccodeupright']=2
factor['toriccode3D']=3

N=args.B


if args.Simple:
    S=np.zeros(N)
    for k in range(3,len(renyi_dataG[0].x),args.Step):
        beta=np.array(renyi_dataG[0].x)[:k]
        y=np.array(renyi_dataG[0].y)[:k]
        av_magn=np.array([float(str(y[i]).split(" ")[0]) for i in range(len(y))])
        av_magn_err=np.array([float(str(y[i]).split(" ")[2]) for i in range(len(y))])

        for i in range(N):
            Bstrap_magn=[np.random.normal(av_magn[j],av_magn_err[j]) for j in range(len(y))]
            S[i]=simps(Bstrap_magn,beta,even='even')
        print beta


        gamma=1./(n-1)*np.array([args.Const-(n*factor[args.System]*renyi_dataG[0].props['L']**args.Dim)*S[l] for l in range(N)])

        if args.outfile==None:
            print beta[-1], np.mean(gamma), np.std(gamma)#, np.mean(SZ2), np.mean(S[0]), np.mean(S[1]), np.mean(S[2])
            
else:
    S=np.zeros((5,N))

    with open(args.outfile,'w') as outfile:
        for k in range(3,len(renyi_dataG[0].x),args.Step):
            beta=np.array(renyi_dataG[0].x)[:k]
            for schemes in renyi_dataG:
                scheme=int(schemes.props['IncNo'])
                y=np.array(schemes.y)[:k]
                av_magn=np.array([float(str(y[i]).split(" ")[0]) for i in range(len(y))])
                av_magn_err=np.array([float(str(y[i]).split(" ")[2]) for i in range(len(y))])

                for i in range(N):
                    Bstrap_magn=[np.random.normal(av_magn[j],av_magn_err[j]+0.001*abs(av_magn[j])) for j in range(len(y))]
                    S[scheme,i]=simps(Bstrap_magn,beta,even='avg')

            if args.alternate:
                gamma=np.array([args.Const-(2*factor[args.System]*schemes.props['L']**args.Dim)*(S[0,l]-S[1,l]-S[2,l]+S[3,l]) for l in range(N)])
            else :
                gamma=1./(n-1)*np.array([args.Const-(n*factor[args.System]*schemes.props['L']**args.Dim)*(S[0,l]-2*S[1,l]+S[2,l]) for l in range(N)])

            if args.outfile==None:
                print beta[-1], np.mean(gamma), np.std(gamma)#, np.mean(SZ2), np.mean(S[0]), np.mean(S[1]), np.mean(S[2])
            else:
                outfile.write("{0:f} {1:f} {2:f} \n".format(beta[-1], np.mean(gamma), np.std(gamma)) )
        outfile.close()
