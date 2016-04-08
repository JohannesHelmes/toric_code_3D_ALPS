import uncertainties as unc
import numpy as np
import pyalps

def resampling(avg,err,N=1000):
    """ Performs resampling of an array of normal distributed data

    Parameters
    ----------
    avg : array_like
            Gaussian distributed average values
    err : array_like
            Errors on values
    N   : int, optional
            Number of samples (default = 1000)

    Returns
    ---------
    samples : ndarray
            resampled values, the result is an array of dimension avg.ndim + 1
    """

    #Needs numpy
    if len(avg) != len(err):
        raise ValueError("Arrays must be of the same dimension.")

    
    return np.array( [ np.random.normal(avg[j],err[j],size=N) if err[j]>0 else avg[j]*np.ones(N) for j in range(len(avg)) ] )


def bootstrap(f,args,argserr,N=1000):
    """ Performs a bootstrapping of a function using normal distributed resampling

    Parameters
    ----------
    f   : callable function
            the function that returns the quantity that should be computed
    args : array_like
            arguments that are passed to the function, i.e. the command f(*args) must be possible
    argserr : array_like
            errors on the argument in the same order as in args. Will be passed to the resampling
    N   : int, optional
            Number of values for resampling(default = 1000)

    Returns
    ---------
    y, yerr : floats
        outcome of the function +/- error
            
    """

    if len(args) != len(argserr):
        raise ValueError("Array of arguments has different size than array of argument errors.")

    sample_array=np.empty((len(args),N))
    for i, (a, aerr) in enumerate(zip(args, argserr)):
        sample_array[i,:]=resampling([a], [aerr],N)

    y_list=np.empty(N)
    for j in xrange(N):
        y_list[j] = f(*sample_array[:,j])


    return np.mean(y_list), np.std(y_list,ddof=1)/np.sqrt(N)

    



def thermodynamic_integration(x,y,yerr,N=1000):
    """ Performs the integration using Simpson's rule after bootstrapping the y values.
    """
    Q=np.zeros(N)

    y_bstrap=bootstrap(y,yerr,N)
    for i in xrange(N):
        Q[i]=simps(y_bstrap[:,i],x,even='avg')

    return Q
    
def read_results_from_file(path, prefix, X, Y, ForEach=None):
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
        dataset= pyalps.loadMeasurements(pyalps.getResultFiles(dirname=path,prefix=prefix),Y)

        if ForEach==None:
            return pyalps.collectXY(dataset, x=X, y=Y)
        else:
            return pyalps.collectXY(dataset, x=X, y=Y, foreach=[ForEach])

    else:
        raise ValueError("Not yet implemented for multiple runs.")

def convert_alps_dataset(dset):
    """ Converts an ALPS dataset to three numpy arrays 

    Parameters
    ----------
    dset       : ALPS dataset

    Returns
    ---------
    X, Y, Yerr : numpy arrays
    """
    X=np.array(dset.x)
    Y = np.array(map(lambda x: float(str(x).split(" ")[0]), dset.y))
    Yerr = np.array(map(lambda x: float(str(x).split(" ")[2]), dset.y))
    return X,Y,Yerr


def convert_alps_dataset_unc(dset):
    """ Converts an ALPS dataset to a numpy array for X, and an uncertainties array for Y +/- Yerr 

    Parameters
    ----------
    dset       : ALPS dataset

    Returns
    ---------
    X       : numpy array
    Y, Yerr : uncertainties array
    """
    X=dset.x
    Y = map(lambda x: unc.ufloat_fromstr(str(x)), dset.y)
    return X,Y
