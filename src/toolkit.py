import numpy as np

def quantile(x, q, weights=None):
    x = np.atleast_1d(x)
    q = np.atleast_1d(q)

    if weights is None:
        return np.percentile(x, 100.0 * q)
    else:
        weights = np.atleast_1d(weights)
        idx = np.argsort(x, axis=0)

        res = []
        for i in range(x.shape[1]):
            sw = weights[idx[:,i]]
            cdf = np.cumsum(sw)[:-1]
            cdf /= cdf[-1]
            cdf = np.append(0, cdf)
            res.append(np.interp(q, cdf, x[idx[:,i],i]))
        return np.array(res).T
        
def get_binned_median(xdata, ydata, bins):
    '''
    return median vals and the standand deviation of the mean for given data and a bin.
    Input:
        xdata: array-like[N,]
        ydata: array-like[N,]
        bins: array-like[Nbin,]
    Output:
        xcs: array-like[Nbin,]
            center of each bin
        medians: array-like[Nbin,]
            median values of ydata in each bin
        stds: array-like[Nbin,]
            stds of the mean of ydata in each bin
    '''
    
    Nbin = len(bins)-1
    medians, stds = [np.zeros(Nbin) for i in range(2)]
    for ibin in range(Nbin):
        idx = (xdata>bins[ibin]) & (xdata<=bins[ibin+1]) & (np.isfinite(xdata)) & (np.isfinite(ydata))
        if np.sum(idx)==0:
            medians[ibin], stds[ibin] = np.nan, np.nan
        else:
            medians[ibin] = np.median(ydata[idx])
            stds[ibin] = np.std(ydata[idx])/np.sqrt(np.sum(idx))
    xcs = (bins[1:]+bins[:-1])/2.
    return xcs, medians, stds