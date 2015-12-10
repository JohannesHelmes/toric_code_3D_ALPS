#!/bin/bash alpspython

import sys,os,time
import pyalps,argparse
import numpy as np
import configparser

config = configparser.ConfigParser()

parser = argparse.ArgumentParser(description='Create input files for ALPS simulations of upright toric code models in 2 and 3 dimensions', epilog='(C) Johannes Helmes 2015')
filegroup=parser.add_argument_group()

filegroup.add_argument('paramfilename')#,type=file)

paramgroup=parser.add_argument_group()


paramgroup.add_argument('--sweeps','-s', type=int,default=10000,help='Number of sweeps, default=10000')
paramgroup.add_argument('--therm','-t', type=int,default=1000,help='Number of thermalization steps, default=1000')
paramgroup.add_argument('--n','-n', type=int,default=2,help='Order of Renyi entropy, default=2')
paramgroup.add_argument('--length','-l', type=int,help='Lattice length',required=True)
paramgroup.add_argument('--lattice','-a', type=int,help='Lattice type, (0 = toric code, 1 = toric code fcr, 2 = toric code 3D), default = 2',default=2)
paramgroup.add_argument('--ExcType','-e', type=int,help='Type of Excitation, (1 = plaquette, 2 = vertex), default = 2',default=2)

paramgroup.add_argument('--beta','-b', nargs='+',type=float,help='beta-1/T, range possible default=2*l')
paramgroup.add_argument('--magnetization','-m', nargs='+',type=float,default=[0.0],help='Magnetization (h), range possible, default=0.0')
paramgroup.add_argument('--infile','-i', help='Prefix of .in.xml')
paramgroup.add_argument('--geofile','-g', type=file, help='File containing the geometry of the increments')
paramgroup.add_argument('--partsize','-p',default=1000, type=int,help='Size subsimulation partition for better distribution on HPC, default=1000')
paramgroup.add_argument('--temper','-T',help='interpret beta as temperature',action='store_true')
paramgroup.add_argument('--verbose','-v',help='Print some extra stuff',action='store_true')
paramgroup.add_argument('--directory','-d',help='Directory of the parameter files, default= /scratch/helmes/simulations/"name of parent directory"/')
args=parser.parse_args()

if args.directory==None:
    dirpath='/scratch/helmes/simulations/'+os.path.split(os.path.split(os.getcwd())[0])[1]
    print dirpath
else:
    dirpath=args.directory

config['DEFAULT'] = { 'sweeps': args.sweeps, 'therm': args.therm, 'n': args.n, 'length': args.length, 'lattice': args.lattice, 'ExcType':args.ExcType, 'beta': args.beta,
        'magnetization': args.magnetization, 'infile': args.infile, 'geofile':args.geofile.name, 'partsize':args.partsize, 'temper':args.temper, 'directory':dirpath} 

with open(args.paramfilename,'w') as configfile:
    timestring=time.strftime("%Y-%m-%d %H:%M:%S \n\n", time.localtime())
    configfile.write("# created on " + timestring)
    config.write(configfile)

beta=args.beta

latticename=["toric code","toric code fcr","toric code 3D"]
dimension=3 if args.lattice==2 else 2

if (beta==None):
    beta=[2*args.length]

if args.geofile==None:
    geo="0"
else:
    geo=os.path.abspath(args.geofile.name)


parms = []
betagrid=beta if len(beta)==3 else np.array((beta[0],beta[0],1))
magngrid=np.array(args.magnetization) if len(args.magnetization)==3 else np.array((args.magnetization[0],args.magnetization[0],1))
if args.verbose:
    print "magngrid ", magngrid

betalist,hlist= np.mgrid[betagrid[0]:betagrid[1]:betagrid[2]*1j,magngrid[0]:magngrid[1]:magngrid[2]*1j]
if args.verbose:
    print "Blist", Blist
if geo=="0":
    Geometry=[""]
else:
    Geometry=np.genfromtxt(geo,dtype=str,skip_header=0)

if Geometry.size==1:
    Geometry=[str(Geometry)]

if len(Geometry[0])!=dimension*args.length**dimension: #first factor is the number of spins of the unit cell
    print "Warning: given length does not fit to geometry file ",len(Geometry[0])," != ",args.length,"^", dimension

if args.verbose:
    print Geometry

for i,IncEl in enumerate(Geometry):

    for beta,h in zip(betalist.flat,hlist.flat):
        if args.verbose:
            print B,beta,h
 
        parms.append(
            {
                'LATTICE_LIBRARY' : "mylatticelib.xml",
                'LATTICE'	: latticename[args.lattice],
                'h'         : float(h),
                'beta'		: 1./float(beta) if args.temper else float(beta),
                'THERMALIZATION': args.therm,
                'SWEEPS'	: args.sweeps,
                'L'		: args.length,
                'n'		    : args.n,
                'ExcType'	: args.ExcType,
                'd'		: dimension,
                'M'		: 100,
                'IncStep' : IncEl,
                'IncNo' : i
     
            }
        )

cuts=np.append(np.arange(0,len(parms),args.partsize),len(parms))
print cuts
#path=os.path.dirname(args.infile)

if not os.path.exists(dirpath):
    os.makedirs(dirpath)
    print dirpath,"created"

for i,[s,e] in enumerate(zip(cuts[:-1],cuts[1:])):
    #fname=path+'/part'+str(i)+'/'+os.path.basename(args.infile)+'.'+str(i)
    fname=dirpath+'/'+os.path.basename(args.infile)+'.'+str(i)
    #os.system('mkdir -p '+path+'/part'+str(i))
    input_file = pyalps.writeInputFiles(fname,parms[s:e])

