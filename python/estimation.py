import constants
import errorbars
import utils
import numpy as np
from scipy.linalg import block_diag
from scipy.signal import butter, filtfilt
from scipy.interpolate import make_lsq_spline


def interpolate(epochs, yaw):
    yaw = filterLowpass(yaw, Wn=0.1, N=8)
    t = []
    k = 10
    t = np.r_[(epochs[0],) * (k + 1),
              t,
              (epochs[-1],) * (k + 1)]
    f = make_lsq_spline(epochs, yaw, t, k)

    return f(epochs)


def filterLowpass(residualSignal, Wn=0.5, N=8):
    """
    This function gets an input signal and filters it using a lowpass Butterworth filter.
    TODO: try multiple filtering types
    :param residualSignal: residual signal to be filtered
    :return: signal after filtering
    """
    b, a = butter(N, Wn, 'low')
    filteredSignal = filtfilt(b, a, residualSignal)
    return filteredSignal


def filterResiduals(residualData):
    """
    This function applies the filtering to data stored in a residual data dictionary.
    :param residualData: dictionary containing residual signals as read from .res file
    :return: dictionary containing filtered residual signals
    """
    stationsList = list(residualData.keys())

    for station_index in stationsList:
        residualData[station_index]['rk'] = filterLowpass(residualData[station_index]['rk'])
    return residualData


def computeBlock1DesignMatrix(residualData, orbitData, stationsData, epoch):
    """
    This function computes the Least squares design matrix for the first  model.
    This matrix has 2 columns and N_stations rows and is computed for only one epoch.
    :param residualData: dictionary containing residual signals as read from .res file
    :param orbitData: dictionary containing data read from .orb file, used to compute azimuth and nadir angles for the design matrix
    :param stationsData: dictionary containing stations data (index, name and position)
    :param epoch: single epoch in MJD format
    :return: A and rk - the LSE design matrix and measurements vector containing residuals read from residualData
    """
    rk = []
    ak = []
    nk = []
    for station_index in residualData.keys():
        if epoch in residualData[station_index]['mjd']:
            epoch_index = residualData[station_index]['mjd'].index(epoch)

            r_sat, v_sat = utils.getSatellitePositionVelocity(orbitData, epoch)
            azi, nad = utils.getAzimuthNadirSatellite(r_sat, v_sat, station_index, stationsData)

            rk.append(residualData[station_index]['rk'][epoch_index])
            ak.append(azi)
            nk.append(nad)
    A = np.zeros((len(rk), 2))
    for i in range(0, len(rk)):
        A[i][0] = np.sin(nk[i]) * np.sin(ak[i])
        A[i][1] = np.sin(nk[i]) * np.cos(ak[i])
    return A, rk


def computeBlock2DesignMatrix(yaw): # TODO comment code for God's sake
    """
    This function computes the restrictions block which is concatenated underneath block 1.
    The block needs a initial estimated yaw value.
    :param yaw: initial yaw value used for restriction
    :return: block 2 for a certain epoch (contains only 2 elements) and measurement to be concatenated
    at the end of the residual vector: A2 and r2
    """
    A2 = np.array([[2 * constants.PCO_x * np.cos(np.deg2rad(yaw)),
                    2 * constants.PCO_x * np.sin(np.deg2rad(yaw))]])
    r2 = np.array([2 * (constants.PCO_x ** 2)])
    return A2, r2


def solveLSE(A, rk):
    """
    This function gets a design matrix and a residual vector (for one epoch) and solves the LSE estimation.
    This is used for model 1, as the other models require the design matrix to be modified.
    :param A: design matrix for LSE system
    :param rk: measurements vector (residual)
    :return: yaw estimated from LSE solution
    """
    AtA = A.transpose().dot(A)
    x = np.linalg.inv(AtA).dot(A.transpose()).dot(rk)
    yaw = np.rad2deg(np.arctan2(x[1], x[0]))

    errors = errorbars.computeFormalError(A, rk, x)
    yaw_error = np.rad2deg(np.arctan2(errors[1], errors[0]))

    return yaw, yaw_error


def solveLSE2(A, rk, y1):
    """
    This function gets a design matrix and a residual vector (for one epoch) and solves the LSE estimation.
    This is used for model 2, as the other models require the design matrix to be modified.
    :param A: design matrix for LSE system
    :param rk: measurements vector (residual)
    :return: yaw estimated from LSE solution
    """
    #y1, _ = solveLSE(A, rk)
    A2 = np.append(A, [[2 * constants.PCO_x * np.cos(np.deg2rad(y1)),
                        2 * constants.PCO_x * np.sin(np.deg2rad(y1))]], axis=0)
    r2 = np.append(rk, [2 * (constants.PCO_x ** 2)], axis=0)
    AtA = A2.transpose().dot(A2)
    x = np.linalg.inv(AtA).dot(A2.transpose()).dot(r2)

    errors = errorbars.computeFormalError(A2, r2, x)
    yaw_error = np.rad2deg(np.arctan2(errors[1], errors[0]))

    y2 = np.rad2deg(np.arctan2(x[1], x[0]))
    return y2, yaw_error


def solveLSE3(A, rk, y1):
    """
    This function gets a design matrix and a residual vector (for one epoch) and solves the LSE estimation.
    This is used for model 3, as the other models require the design matrix to be modified.
    :param A: design matrix for LSE system
    :param rk: measurements vector (residual)
    :return: yaw estimated from LSE solution
    """
    #y1, _ = solveLSE(A, rk)
    A2 = np.append(A, [[2 * constants.PCO_x * np.cos(np.deg2rad(y1)),
                        2 * constants.PCO_x * np.sin(np.deg2rad(y1))]], axis=0)
    b = np.ones((len(rk), 1))
    b = np.concatenate((b, np.zeros((1, 1))), axis=0)

    A2 = np.concatenate((A2, b), axis=1)

    r2 = np.append(rk, [2 * (constants.PCO_x ** 2)], axis=0)
    AtA = A2.transpose().dot(A2)
    x = np.linalg.inv(AtA).dot(A2.transpose()).dot(r2)

    errors = errorbars.computeFormalError(A2, r2, x)
    yaw_error = np.rad2deg(np.arctan2(errors[1], errors[0]))

    y2 = np.rad2deg(np.arctan2(x[1], x[0]))
    return y2, yaw_error


def solveLSE3ClockCorrections(A, rk, y1):
    A2 = np.append(A, [[2 * constants.PCO_x * np.cos(np.deg2rad(y1)),
                        2 * constants.PCO_x * np.sin(np.deg2rad(y1))]], axis=0)
    b = np.ones((len(rk), 1))
    b = np.concatenate((b, np.zeros((1, 1))), axis=0)

    A2 = np.concatenate((A2, b), axis=1)

    r2 = np.append(rk, [2 * (constants.PCO_x ** 2)], axis=0)
    AtA = A2.transpose().dot(A2)
    x = np.linalg.inv(AtA).dot(A2.transpose()).dot(r2)

    errors = errorbars.computeFormalError(A2, r2, x)
    clock_error = errors[2]

    return x[2], clock_error


def solveLSEModel1(residualData, orbitData, stationsData):
    """
    This function implements a simple LSE model on residual data to estimate yaw angles.
    :param residualData: dictionary containing data read from .res file
    :param orbitData: dictionary containing data read from .orb file, used to compute azimuth and nadir angles for the design matrix
    :param stationsData: dictionary containing stations data (index, name and position)
    :return: estimated yaw angles in degrees and their corresponding epochs
    """
    epochs = utils.getEpochsArray(residualData)
    yaw = []
    errors = []

    for epoch in epochs:
        A, rk = computeBlock1DesignMatrix(residualData, orbitData, stationsData, epoch)
        y, error = solveLSE(A, rk)

        yaw.append(y)
        errors.append(error)
    yaw = shiftYawPlacement(yaw, orbitData['yaw'])
    return (epochs, yaw, errors)


def solveLSEModel2(residualData, orbitData, stationsData):
    """
    This function implements a simple LSE model on residual data to estimate yaw angles, while using a restraint condition.
    :param residualData: dictionary containing data read from .res file
    :param orbitData: dictionary containing data read from .orb file, used to compute azimuth and nadir angles for the design matrix
    :param stationsData: dictionary containing stations data (index, name and position)
    :return: estimated yaw angles in degrees and their corresponding epochs
    """
    epochs = utils.getEpochsArray(residualData) #TODO add initial conditions outside of solveLSE2 - to implement filtering !!!
    yaw = []
    errors = []
    _, yaw_init, _ = solveLSEModel1(residualData, orbitData, stationsData)
    yaw_init = shiftYawPlacement(yaw_init, orbitData['yaw'])
    #yaw_init = interpolate(epochs, yaw_init)

    for epoch in epochs:
        A, rk = computeBlock1DesignMatrix(residualData, orbitData, stationsData, epoch)
        y2, error = solveLSE2(A, rk, yaw_init[epochs.index(epoch)])

        yaw.append(y2)
        errors.append(error)

    yaw = shiftYawPlacement(yaw, orbitData['yaw'])
    return (epochs, yaw, errors)


def solveLSEModel3(residualData, orbitData, stationsData):
    """
    This function implements a simple LSE model on residual data to estimate yaw angles, while using a restraint condition
    and also estimating satellite clock corrections. TODO add satellite clocks to output and plots.
    :param residualData: dictionary containing data read from .res file
    :param orbitData: dictionary containing data read from .orb file, used to compute azimuth and nadir angles for the design matrix
    :param stationsData: dictionary containing stations data (index, name and position)
    :return: estimated yaw angles in degrees and their corresponding epochs
    """
    epochs = utils.getEpochsArray(residualData)
    epochs_final = []
    yaw = []
    errors = []

    _, yaw_init, _ = solveLSEModel2(residualData, orbitData, stationsData)
    yaw_init = shiftYawPlacement(yaw_init, orbitData['yaw'])
    # yaw_init = interpolate(epochs, yaw_init)

    for epoch in epochs:
        A, rk = computeBlock1DesignMatrix(residualData, orbitData, stationsData, epoch)
        try:
            y2, error = solveLSE3(A, rk, yaw_init[epochs.index(epoch)])

            yaw.append(y2)
            epochs_final.append(epoch)
            errors.append(error)

        except np.linalg.LinAlgError:
            print("Singular matrix at epoch {}.".format(epoch))
    yaw = shiftYawPlacement(yaw, orbitData['yaw'])
    return (epochs_final, yaw, errors)


def getCLockCorrections(residualData, orbitData, stationsData):
    epochs = utils.getEpochsArray(residualData)
    epochs_final = []
    times_corr = []
    errors = []

    _, yaw_init, _ = solveLSEModel2(residualData, orbitData, stationsData)
    yaw_init = shiftYawPlacement(yaw_init, orbitData['yaw'])
    # yaw_init = interpolate(epochs, yaw_init)

    for epoch in epochs:
        A, rk = computeBlock1DesignMatrix(residualData, orbitData, stationsData, epoch)
        try:
            tc, error = solveLSE3ClockCorrections(A, rk, yaw_init[epochs.index(epoch)])

            times_corr.append(tc)
            epochs_final.append(epoch)
            errors.append(error)

        except np.linalg.LinAlgError:
            print("Singular matrix at epoch {}.".format(epoch))
    return (epochs_final, times_corr, errors)


def computeBlock1DesignMatrixWindow(residualData, orbitData, stationsData, startEpoch, endEpoch):
    """
    This function computes the first block of the grand design matrix used in model 4. This model estimates yaw angles
    over multiple epochs and uses multiple stations - therefore the first block which is computed here
    will have smaller blocks on the main diagonal (identical to those used for the first model).
    The final block 1 returned from this funciton will have 2*N_epochs columns and N_epochs*N_stations rows.
    :param residualData: dictionary containing data read from .res file
    :param orbitData: dictionary containing data read from .orb file, used to compute azimuth and nadir angles for the design matrix
    :param stationsData: dictionary containing stations data (index, name and position)
    :param startEpoch: float variable containing the first epoch, in MJD format, for the interval over which the
    block matrix is computed
    :param endEpoch: float variable containing the last epoch, in MJD format, for the interval over which the
    block matrix is computed
    :return: Block 1 matrix and residual data
    """
    trimmedResData = utils.trimResDataWindow(residualData, startEpoch, endEpoch)
    epochs = utils.getEpochsArray(trimmedResData)

    grandDesignMatrixBlock1 = np.empty((0, 0), float)
    Rk = np.empty((0, 0), float)
    for epoch in epochs:
        A, rk = computeBlock1DesignMatrix(trimmedResData, orbitData, stationsData, epoch)
        grandDesignMatrixBlock1 = block_diag(grandDesignMatrixBlock1, A)
        try:
            Rk = np.concatenate((Rk, np.array(rk)), axis=0)
        except:
            Rk = np.array(rk)
    return (grandDesignMatrixBlock1, Rk)


def computeBlock2DesignMatrixWindow(residualData, orbitData, stationsData, startEpoch, endEpoch, yaw):
    """
    This function computes the restriction matrices for multiple epochs (contained in the interval between
    startEpoch and endEpoch) and concatenates them in a block-diagonal matrix. The final block has N_epochs
    rows and 2*N_epochs columns.
    """
    trimmedResData = utils.trimResDataWindow(residualData, startEpoch, endEpoch)
    epochs = utils.getEpochsArray(trimmedResData)

    #epochs, _ = solveLSEModel3(trimmedResData, orbitData, stationsData)
    yaw_init = yaw[orbitData['mjd'].index(startEpoch): orbitData['mjd'].index(endEpoch)]
    #yaw_init = filterLowpass(yaw_init, Wn=0.1, N=4)

    grandDesignMatrixBlock2 = np.empty((0, 0), float)
    Rk = np.empty((0, 0), float)
    for epoch in epochs:

        A, rk = computeBlock2DesignMatrix(yaw_init[epochs.index(epoch)])
        grandDesignMatrixBlock2 = block_diag(grandDesignMatrixBlock2, A)
        try:
            Rk = np.concatenate((Rk, np.array(rk)), axis=0)
        except:
            Rk = np.array(rk)
    return (grandDesignMatrixBlock2, Rk)


def computeBlock3DesignMatrixWindow(residualData, startEpoch, endEpoch):
    """
    This function computes the 3rd block for the fourth model LSE estimator. This block has the values
    (coefficients) multiplied with satellite clock errors and stations biases.
    """
    trimmedResData = utils.trimResDataWindow(residualData, startEpoch, endEpoch)
    epochs = utils.getEpochsArray(trimmedResData)
    nr_st = len(trimmedResData.keys())

    grandDesignMatrixBlock3 = np.empty((0, 0), float)
    for epoch in epochs:
        C = np.zeros((nr_st, len(epochs)))
        C[:, epochs.index(epoch)] = 1
        C = np.concatenate((C, 0.1*np.eye(nr_st)), axis=1)
        try:
            grandDesignMatrixBlock3 = np.concatenate((grandDesignMatrixBlock3, C), axis=0)
        except:
            grandDesignMatrixBlock3 = C
    return grandDesignMatrixBlock3


def computeBlock4DesignMatrixWindow(residualData, startEpoch, endEpoch):
    """
    Computes the final block to be concatenated to the grand design matrix for the fourth LSE model.
    This block contains only zero values over a dimension of N_epochs*(N_epochs*N_stations).
    """
    trimmedResData = utils.trimResDataWindow(residualData, startEpoch, endEpoch)
    epochs = utils.getEpochsArray(trimmedResData)
    nr_st = len(trimmedResData.keys())

    grandDesignMatrixBlock4 = np.zeros((len(epochs), len(epochs)+nr_st))
    return grandDesignMatrixBlock4


def computeGrandDesignMatrix(residualData, orbitData, stationsData, startEpoch, endEpoch, yaw_init):
    """
    Computes the grand design matrix using the four smaller blocks from the previous functions.
    This matrix corresponds to the fourth estimation model, which takes into account integer ambiguities,
    satellite clock corrections and multiple epochs.
    """
    A1, rk = computeBlock1DesignMatrixWindow(residualData, orbitData, stationsData, startEpoch, endEpoch)
    A2, rk2 = computeBlock2DesignMatrixWindow(residualData, orbitData, stationsData, startEpoch, endEpoch, yaw_init)
    A3 = computeBlock3DesignMatrixWindow(residualData, startEpoch, endEpoch)
    A4 = computeBlock4DesignMatrixWindow(residualData, startEpoch, endEpoch)

    grandDesignMatrix = np.concatenate((A1, A2), axis=0)
    B = np.concatenate((A3, A4), axis=0)
    grandDesignMatrix = np.concatenate((grandDesignMatrix, B), axis=1)
    Rk = np.concatenate((rk, rk2), axis=0)

    #Ne = 360
    #grandDesignMatrix = np.concatenate((A1, A2), axis=0)
    #B = np.concatenate((A3[:, 0:Ne], A4[:, 0:Ne]), axis=0)
    #grandDesignMatrix = np.concatenate((grandDesignMatrix, B), axis=1)
    #Rk = np.concatenate((rk, rk2), axis=0)

    return grandDesignMatrix, Rk


def solveLSE4(residualData, orbitData, stationsData, startEpoch, endEpoch, yaw_init):
    """
    This function solves the fourth model and outputs an array of yaw values alongside their corresponding epochs.
    All yaw values correspond to the interval between startEpoch and endEpoch.
    """
    A, rk = computeGrandDesignMatrix(residualData, orbitData, stationsData, startEpoch, endEpoch, yaw_init)

    trimmedResData = utils.trimResDataWindow(residualData, startEpoch, endEpoch)
    epochs = utils.getEpochsArray(trimmedResData)
    Ne = len(epochs)

    AtA = A.transpose().dot(A)
    x = np.linalg.inv(AtA).dot(A.transpose()).dot(rk)

    yaw = []
    yaw_error = []

    errors = errorbars.computeFormalError(A, rk, x)
    for i in range(0, Ne*2, 2):
        x0 = x[i]
        y0 = x[i+1]
        yaw.append(np.rad2deg(np.arctan2(y0, x0)))
        yaw_error.append(np.rad2deg(np.arctan2(errors[i+1], errors[i])))

    return epochs, yaw, yaw_error


def solveLSEModel4(residualData, orbitData, stationsData, Ne=5):
    """
    This function iterates through groups of epochs (grouped in bins with Ne epochs each).
    For each bin the fourth model estimation is applied. This can also be applied for all available epochs
    at once yielding more precise results than applying it for smaller time intervals.
    """
    #epochs, yaw_init, _ = solveLSEModel3(residualData, orbitData, stationsData)
    #yaw_init = filterLowpass(yaw_init, Wn=0.1, N=4)

    epochs, yaw_init, _ = solveLSEModel3(residualData, orbitData, stationsData)
    yaw_init = shiftYawPlacement(yaw_init, orbitData['yaw'])
    #yaw_init = interpolate(epochs, yaw_init)

    yaw = []
    epoch_final = []
    errors = []
    for epoch in epochs[::Ne]:
        try:
            startEpoch = epoch
            endEpoch = epochs[epochs.index(epoch)+Ne]
            e, y, error = solveLSE4(residualData, orbitData, stationsData, startEpoch, endEpoch, yaw_init)

            yaw = yaw + y
            epoch_final = epoch_final + e
            errors = errors + error
        except (np.linalg.LinAlgError, IndexError):
            print("Singular A matrix for epoch {}".format(epoch))
    yaw = shiftYawPlacement(yaw, orbitData['yaw'])
    return epoch_final, yaw, errors


def shiftYawPlacement(yaw_estimated, yaw_nominal, Ne=3): # TODO shifting function AND polyfit for initial yaw values
    for yaw in yaw_estimated:
        index = yaw_estimated.index(yaw)
        start_value = yaw_nominal[yaw_estimated.index(yaw)]
        if np.sqrt((yaw - start_value) ** 2) > 200 and index <= Ne:
            yaw_estimated[index] = yaw - 360 * np.sign(yaw_estimated[index])
        else:
            median_value = np.median(yaw_estimated[index-Ne:index])
            for _ in range(0, 2):
                if np.sqrt((yaw - median_value) ** 2) > 200:
                    yaw_estimated[index] = yaw + 360 * np.sign(yaw_estimated[index-1])
                    #yaw_estimated[index] = yaw + 360 * np.sign(median_value)
                    #if np.sqrt((yaw - median_value) ** 2) > 200 or np.sqrt((yaw - yaw_estimated[index - 1]) ** 2) > 200:
                    #    yaw_estimated[index] = yaw - 2 * 360 * np.sign(yaw_estimated[index - 1])
    return yaw_estimated



