import constants
import utils
import numpy as np
from scipy.linalg import block_diag
from scipy.signal import butter, filtfilt


def filterLowpass(residualSignal):
    """
    This function gets an input signal and filters it using a lowpass Butterworth filter.
    TODO: try multiple filtering types
    :param residualSignal: residual signal to be filtered
    :return: signal after filtering
    """
    Wn = 0.4
    N = 5
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


def solveLSE(A, rk):
    AtA = A.transpose().dot(A)
    x = np.linalg.inv(AtA).dot(A.transpose()).dot(rk)
    yaw = -np.rad2deg(np.arctan2(x[1], x[0]))
    return yaw


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
        y1 = solveLSE(A, rk)
        # y1 = orbitData['yaw'][orbitData['mjd'].index(epoch)] # TODO see if using nominal yaw angle as starting point yields better results

        AtA = A.transpose().dot(A)
        x = np.linalg.inv(AtA).dot(A.transpose()).dot(rk)
        y1 = np.rad2deg(np.arctan2(x[1], x[0]))

        A2 = np.append(A, [[2 * constants.PCO_x * np.cos(np.deg2rad(y1)),
                            2 * constants.PCO_x * np.sin(np.deg2rad(y1))]], axis=0)
        r2 = np.append(rk, [2 * (constants.PCO_x ** 2)], axis=0)
        AtA = A2.transpose().dot(A2)
        x = np.linalg.inv(AtA).dot(A2.transpose()).dot(r2)

        y2 = np.rad2deg(np.arctan2(x[1], x[0]))
        yaw.append(-y2) # TODO see why the algorithm seems to estimate -yaw instead of yaw ?!

    return (epochs, yaw)


def solveLSEModel3(residualData, orbitData, stationsData):
    """
    This function implements a simple LSE model on residual data to estimate yaw angles, while using a restraint condition.
    :param residualData: dictionary containing data read from .res file
    :param orbitData: dictionary containing data read from .orb file, used to compute azimuth and nadir angles for the design matrix
    :param stationsData: dictionary containing stations data (index, name and position)
    :return: estimated yaw angles in degrees and their corresponding epochs
    """
    epochs = []
    for station_index in residualData.keys():
        epochs = epochs + list(set(residualData[station_index]['mjd']) - set(epochs))
    epochs.sort()
    yaw = []

    for epoch in epochs:
        A, rk = computeBlock1DesignMatrix(residualData, orbitData, stationsData, epoch)
        y1 = solveLSE(A, rk)
        # y1 = orbitData['yaw'][orbitData['mjd'].index(epoch)] # TODO see if using nominal yaw angle as starting point yields better results

        A2 = np.append(A, [[2 * constants.PCO_x * np.cos(y1), 2 * constants.PCO_x * np.sin(y1)]], axis=0)

        b = np.ones((len(rk), 1))
        b = np.concatenate((b, np.zeros((1, 1))), axis=0)

        A2 = np.concatenate((A2, b), axis=1)
        r2 = np.append(rk, [2*((constants.PCO_x)**2)], axis=0)
        y2 = solveLSE(A2, r2)

        yaw.append(y2) # TODO see why the algorithm seems to estimate -yaw instead of yaw ?!

    return (epochs, yaw)


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
    (epochs_init, yaw_init) = solveLSEModel2(trimmedResData, orbitData, stationsData)

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
        B = np.array([2*constants.PCO_x*np.cos(yaw_init[epochs_init.index(epoch)]),
                      2*constants.PCO_x*np.sin(yaw_init[epochs_init.index(epoch)])])
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



def solveLSEModel4WithRestraint(residualData, orbitData, stationsData, startEpoch, endEpoch):
    # trimmedResData = utils.trimResDataWindow(residualData, startEpoch, endEpoch)
    # (epochs, _) = solveLSEModel1WithRestraint(trimmedResData, orbitData, stationsData)

    (A, rk, epochs) = computeGrandDesignMatrixWindow(residualData, orbitData, stationsData, startEpoch, endEpoch)

    AtA = A.transpose().dot(A)
    X = np.linalg.inv(AtA).dot(A.transpose()).dot(rk)
    yaw = []

    for epoch in epochs:
        i = epochs.index(epoch)
        x = X[2*i: (2*i)+2]
        y = np.rad2deg(np.arctan2(x[1], x[0]))
        yaw.append(y)

    return (epochs, yaw)



def solveLSEModel4WithRestraintWindow(residualData, orbitData, stationsData, Ne=10):
    epochs = utils.getEpochsArray(residualData)
    epochs_final = []
    yaw = []
    for epoch in epochs[1:-Ne-1]:
        index = epochs.index(epoch)
        startEpoch = epoch

        endEpoch = epochs[index+Ne]
        try:
            e, y = solveLSEModel4WithRestraint(residualData, orbitData, stationsData, startEpoch, endEpoch)

            epochs_final.append(e[5])
            yaw.append(y[5])

            print(epochs_final)
            print(yaw)
        except:
            pass

    return (epochs_final, yaw)