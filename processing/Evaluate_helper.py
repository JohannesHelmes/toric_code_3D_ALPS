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

    y_bstrap=bootstrap(y,yerr,N)
    for i in xrange(N):
        Q[i]=simps(y_bstrap[:,i],x,even='avg')

    return Q
    

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
