import constants
import numpy as np
import utils


def computeFormalError(A, r, x):
    """
    This function computes the formal errors (error bars) for a LSE estimation linear system.
    It uses as inputs the design matrix A, the data vector r (containing residuals)
    and the estimated variables vector x.
    """
    n_p = len(x) # number of unknowns - number of estimated variables
    n = len(r)  # number of observations (data samples) - per epoch
    e = r - A.dot(x) # estimation residuals vector

    m = (e.dot(e)) / (n - n_p)
    m0 = np.sqrt(m * np.sign(m))
    Q = np.linalg.inv(A.transpose().dot(A))

    eps = []
    diag = list(np.diagonal(Q))
    for item in diag:
        e = m0 * np.sqrt(item*np.sign(item))
        eps.append(e)

    return eps
