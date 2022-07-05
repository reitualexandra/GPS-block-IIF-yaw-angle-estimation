import constants
import utils
import numpy as np
from scipy.linalg import block_diag
from scipy.signal import butter, filtfilt






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
    rk = []
    ak = []
    nk = []
    for station_index in residualData.keys():
        if epoch in residualData[station_index]['mjd']:
            epoch_index = residualData[station_index]['mjd'].index(epoch)

            r_sat, v_sat = utils.getSatellitePositionVelocity(orbitData, epoch)
            azi, nad = utils.getAzimuthNadirSatellite(r_sat, v_sat, station_index, stationsData)
            # azi, ele = utils.getAzimuthElevationTopocentric(r_sat, station_index, stationsData)

            rk.append(residualData[station_index]['rk'][epoch_index])
            ak.append(azi)
            nk.append(nad)
    A = np.zeros((len(rk), 2))
    for i in range(0, len(rk)):
        A[i][0] = np.sin(nk[i]) * np.sin(ak[i])
        A[i][1] = np.sin(nk[i]) * np.cos(ak[i])
    return A, rk


def computeBlock2DesignMatrix(residualData, orbitData, stationsData, epoch): # TODO comment code for God's sake
    A, rk = computeBlock1DesignMatrix(residualData, orbitData, stationsData, epoch)
    y1 = solveLSE3(A, rk)

    A2 = np.array([[2 * constants.PCO_x * np.cos(np.deg2rad(y1)),
                    2 * constants.PCO_x * np.sin(np.deg2rad(y1))]])
    r2 = np.array([2 * (constants.PCO_x ** 2)])
    return A2, r2


def solveLSE(A, rk):
    AtA = A.transpose().dot(A)
    x = np.linalg.inv(AtA).dot(A.transpose()).dot(rk)
    yaw = np.rad2deg(np.arctan2(x[1], x[0]))
    return yaw


def solveLSE2(A, rk):
    y1 = solveLSE(A, rk)
    A2 = np.append(A, [[2 * constants.PCO_x * np.cos(np.deg2rad(y1)),
                        2 * constants.PCO_x * np.sin(np.deg2rad(y1))]], axis=0)
    r2 = np.append(rk, [2 * (constants.PCO_x ** 2)], axis=0)
    AtA = A2.transpose().dot(A2)
    x = np.linalg.inv(AtA).dot(A2.transpose()).dot(r2)

    y2 = np.rad2deg(np.arctan2(x[1], x[0]))
    return y2


def solveLSE3(A, rk):
    y1 = solveLSE(A, rk)
    A2 = np.append(A, [[2 * constants.PCO_x * np.cos(np.deg2rad(y1)),
                        2 * constants.PCO_x * np.sin(np.deg2rad(y1))]], axis=0)
    b = np.ones((len(rk), 1))
    b = np.concatenate((b, np.zeros((1, 1))), axis=0)

    A2 = np.concatenate((A2, b), axis=1)

    r2 = np.append(rk, [2 * (constants.PCO_x ** 2)], axis=0)
    AtA = A2.transpose().dot(A2)
    x = np.linalg.inv(AtA).dot(A2.transpose()).dot(r2)

    y2 = np.rad2deg(np.arctan2(x[1], x[0]))
    return y2


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

    for epoch in epochs:
        A, rk = computeBlock1DesignMatrix(residualData, orbitData, stationsData, epoch)
        y = solveLSE(A, rk)

        try:
            if np.sqrt((y-yaw[-1])**2) > 250:
                y = y - 360*np.sign(y)
        except IndexError:
            pass

        yaw.append(y) # TODO see why the algorithm seems to estimate -yaw instead of yaw ?!
        # yaw.append(np.rad2deg(np.arctan(x[1]/ x[0])))
    return (epochs, yaw)


def solveLSEModel2(residualData, orbitData, stationsData):
    """
    This function implements a simple LSE model on residual data to estimate yaw angles, while using a restraint condition.
    :param residualData: dictionary containing data read from .res file
    :param orbitData: dictionary containing data read from .orb file, used to compute azimuth and nadir angles for the design matrix
    :param stationsData: dictionary containing stations data (index, name and position)
    :return: estimated yaw angles in degrees and their corresponding epochs
    """
    epochs = utils.getEpochsArray(residualData)
    yaw = []

    for epoch in epochs:
        A, rk = computeBlock1DesignMatrix(residualData, orbitData, stationsData, epoch)

        y2 = solveLSE2(A, rk)

        try:
            if np.sqrt((y2 - yaw[-1]) ** 2) > 250:
                y2 = y2 - 360 * np.sign(y2)
        except IndexError:
            pass

        yaw.append(y2)

    return (epochs, yaw)


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
    yaw = []

    for epoch in epochs:
        A, rk = computeBlock1DesignMatrix(residualData, orbitData, stationsData, epoch)
        y2 = solveLSE3(A, rk)

        try:
            if np.sqrt((y2 - yaw[-1]) ** 2) > 250:
                y2 = y2 - 360 * np.sign(y2)
        except IndexError:
            pass

        yaw.append(y2)

    return (epochs, yaw)


def computeBlock1DesignMatrixWindow(residualData, orbitData, stationsData, startEpoch, endEpoch):
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


def computeBlock2DesignMatrixWindow(residualData, orbitData, stationsData, startEpoch, endEpoch):
    trimmedResData = utils.trimResDataWindow(residualData, startEpoch, endEpoch)
    epochs = utils.getEpochsArray(trimmedResData)

    grandDesignMatrixBlock2 = np.empty((0, 0), float)
    Rk = np.empty((0, 0), float)
    for epoch in epochs:
        A, rk = computeBlock2DesignMatrix(trimmedResData, orbitData, stationsData, epoch)
        grandDesignMatrixBlock2 = block_diag(grandDesignMatrixBlock2, A)
        try:
            Rk = np.concatenate((Rk, np.array(rk)), axis=0)
        except:
            Rk = np.array(rk)
    return (grandDesignMatrixBlock2, Rk)


def computeBlock3DesignMatrixWindow(residualData, startEpoch, endEpoch):
    trimmedResData = utils.trimResDataWindow(residualData, startEpoch, endEpoch)
    epochs = utils.getEpochsArray(trimmedResData)
    nr_st = len(trimmedResData.keys())

    grandDesignMatrixBlock3 = np.empty((0, 0), float)
    for epoch in epochs:
        C = np.zeros((nr_st, len(epochs)))
        C[:, epochs.index(epoch)] = 1
        C = np.concatenate((C, np.eye(nr_st)), axis=1)
        try:
            grandDesignMatrixBlock3 = np.concatenate((grandDesignMatrixBlock3, C), axis=0)
        except:
            grandDesignMatrixBlock3 = C
    return grandDesignMatrixBlock3


def computeBlock4DesignMatrixWindow(residualData, startEpoch, endEpoch):
    trimmedResData = utils.trimResDataWindow(residualData, startEpoch, endEpoch)
    epochs = utils.getEpochsArray(trimmedResData)
    nr_st = len(trimmedResData.keys())

    grandDesignMatrixBlock4 = np.zeros((len(epochs), len(epochs)+nr_st))
    return grandDesignMatrixBlock4


def computeGrandDesignMatrix(residualData, orbitData, stationsData, startEpoch, endEpoch):
    A1, rk = computeBlock1DesignMatrixWindow(residualData, orbitData, stationsData, startEpoch, endEpoch)
    A2, rk2 = computeBlock2DesignMatrixWindow(residualData, orbitData, stationsData, startEpoch, endEpoch)
    A3 = computeBlock3DesignMatrixWindow(residualData, startEpoch, endEpoch)
    A4 = computeBlock4DesignMatrixWindow(residualData, startEpoch, endEpoch)

    grandDesignMatrix = np.concatenate((A1, A2), axis=0)
    B = np.concatenate((A3, A4), axis=0)
    grandDesignMatrix = np.concatenate((grandDesignMatrix, B), axis=1)
    Rk = np.concatenate((rk, rk2), axis=0)

    return grandDesignMatrix, Rk


def solveLSE4(residualData, orbitData, stationsData, startEpoch, endEpoch):
    A, rk = computeGrandDesignMatrix(residualData, orbitData, stationsData, startEpoch, endEpoch)

    trimmedResData = utils.trimResDataWindow(residualData, startEpoch, endEpoch)
    epochs = utils.getEpochsArray(trimmedResData)
    Ne = len(epochs)

    AtA = A.transpose().dot(A)
    x = np.linalg.inv(AtA).dot(A.transpose()).dot(rk)

    yaw = []
    for i in range(0, Ne*2, 2):
        x0 = x[i]
        y0 = x[i+1]
        yaw.append(np.rad2deg(np.arctan2(y0, x0)))

    return epochs, yaw


def final4(residualData, orbitData, stationsData, Ne=5):
    epochs = utils.getEpochsArray(residualData)

    yaw = []
    epoch_final = []
    for epoch in epochs[:-Ne:Ne]:
        try:
            startEpoch = epoch
            endEpoch = epochs[epochs.index(epoch)+Ne]
            e, y = solveLSE4(residualData, orbitData, stationsData, startEpoch, endEpoch)

            yaw.append(y)
            epoch_final.append(e)
        except np.linalg.LinAlgError:
            print("Singular A matrix for epoch {}".format(epoch))

    return epoch_final, yaw


def computeGrandDesignMatrixWindow(residualData, orbitData, stationsData, startEpoch, endEpoch):
    """
    This function computes the grand design matrix for N_stations stations and N_epochs epochs.
    :param residualData: dictionary containing residual data read from .res file
    :param orbitData: dictionary containing orbit data read from .orb file
    :param stationsData: dictionary containing stations data (index, name and position)
    :param startEpoch: epoch at which estimation starts
    :param endEpoch: epoch at which estimation ends
    :return: grand design matrix and data vector with added constraints (fake observations)
    """

    # Get residual data contained between startEpoch and endEpoch into a dictionary format
    trimmedResData = utils.trimResDataWindow(residualData, startEpoch, endEpoch)
    # Get initial epochs and yaw values between startEpoch and endEpoch by using the trimmed residual data
    (epochs_init, yaw_init) = solveLSEModel1(trimmedResData, orbitData, stationsData)
    #yaw_init = filterLowpass(yaw_init, Wn=0.1, N=4)

    # Initialize residual vector which will also contain pseudo observations for restraints
    Rk = np.empty((0, 0), float)

    # Initialize blocks that will make up the grand design matrix
    grandDesignMatrixBlock1 = np.empty((0, 0), float)
    grandDesignMatrixBlock2 = np.empty((0, 0), float)
    grandDesignMatrixBlock3 = np.empty((0, 0), float)


    # For each epoch between startEpoch and EndEpoch get azimuth and nadir. Compute and concatenate buildings blocks.
    for epoch in epochs_init:
        rk = []
        ak = []
        nk = []
        for station_index in trimmedResData.keys():
            # For each epoch construct a list of residual, azimuth and nadir values and concatenate them into arrays
            epoch_index = residualData[station_index]['mjd'].index(epoch)

            r_sat, v_sat = utils.getSatellitePositionVelocity(orbitData, epoch)
            azi, nad = utils.getAzimuthNadirSatellite(r_sat, v_sat, station_index, stationsData)

            rk.append(residualData[station_index]['rk'][epoch_index])
            ak.append(azi)
            nk.append(nad)

        # For each epoch construct a small block with N_stations rows and 2 columns
        A = np.zeros((len(rk), 2))
        for i in range(0, len(rk)):
            A[i][0] = np.sin(np.deg2rad(nk[i])) * np.sin(np.deg2rad(ak[i]))
            A[i][1] = np.sin(np.deg2rad(nk[i])) * np.cos(np.deg2rad(ak[i]))

        # After all small blocks are populated with values put them into a diagonal matrix - the first block
        grandDesignMatrixBlock1 = block_diag(grandDesignMatrixBlock1, A)

        # For each epoch construct the constraints block and concatenate them diagonally
        b1 = 2*constants.PCO_x*np.cos(yaw_init[epochs_init.index(epoch)])
        b2 = 2*constants.PCO_x*np.sin(yaw_init[epochs_init.index(epoch)])
        B = np.array([b1, b2])
        grandDesignMatrixBlock2 = block_diag(grandDesignMatrixBlock2, B)

        # Create the right side block with ones on columns corresponding to epoch + eye matrix on the right
        C = np.zeros((len(rk), len(epochs_init)))
        C[:, epochs_init.index(epoch)] = 1
        C = np.concatenate((C, np.eye(len(rk))), axis=1)
        try:
            grandDesignMatrixBlock3 = np.concatenate((grandDesignMatrixBlock3, C), axis=0)
            Rk = np.concatenate((Rk, np.array(rk)), axis=0)
        except:
            grandDesignMatrixBlock3 = C
            Rk = np.array(rk)

    # Add constraints to residual data vector
    Rk = np.concatenate((Rk, 2*(constants.PCO_x**2)*np.ones(len(epochs_init))), axis=0)
    grandDesignMatrix = np.concatenate((grandDesignMatrixBlock1, grandDesignMatrixBlock2), axis=0)
    grandDesignMatrixBlock3 = np.concatenate((grandDesignMatrixBlock3, np.zeros((len(epochs_init), len(epochs_init)+len(rk)))), axis=0)
    grandDesignMatrix = np.concatenate((grandDesignMatrix, grandDesignMatrixBlock3), axis=1)


    return (grandDesignMatrix, Rk, epochs_init)



def solveLSEModel4(A, rk, epochs):
    # trimmedResData = utils.trimResDataWindow(residualData, startEpoch, endEpoch)
    # (epochs, _) = solveLSEModel1WithRestraint(trimmedResData, orbitData, stationsData)
    # TODO see why on Earth the fourth estimator is not convergent ?!!! For the third the problem was on LSE solving - implementing a separate function fixed it
    AtA = A.transpose().dot(A)
    X = np.linalg.inv(AtA).dot(A.transpose()).dot(rk)
    yaw = []

    for epoch in epochs:
        i = epochs.index(epoch)
        x = X[2*i: (2*i)+2]
        y = np.rad2deg(np.arctan2(x[1], x[0]))
        yaw.append(y)

    #print("Epochs: {}".format(epochs))
    #print("Yaw: {}\n".format(yaw))

    return epochs, yaw



def solveLSEModel4Window(residualData, orbitData, stationsData, Ne=1):
    epochs = utils.getEpochsArray(residualData)

    epochs_final = []
    yaw = []
    for epoch in epochs[0: - Ne - 1]:
        i_start = epochs.index(epoch)
        startEpoch = epoch

        endEpoch = epochs[i_start + Ne]
        (A, rk, e) = computeGrandDesignMatrixWindow(residualData, orbitData, stationsData,
                                                             startEpoch, endEpoch)
        try:
            _, y = solveLSEModel4(A, rk, e)

            epochs_final.append(e)
            yaw.append(y)
        except np.linalg.LinAlgError:
            print("Singular matrix at epoch {}.".format(epoch))
    return (epochs_final, yaw)