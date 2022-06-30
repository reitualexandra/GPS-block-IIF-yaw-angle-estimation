import constants
import utils
import numpy as np
from scipy.linalg import block_diag


def solveLSEModel1(residualData):
    """

    :param residualData:
    :return:
    """
    epochs = utils.getEpochsArray(residualData)
    yaw = []

    for epoch in epochs:
        rk = []
        ak = []
        nk = []
        for station_index in residualData.keys():
            if epoch in residualData[station_index]['mjd']:
                epoch_index = residualData[station_index]['mjd'].index(epoch)

                rk.append(residualData[station_index]['rk'][epoch_index])
                ak.append(residualData[station_index]['azi'][epoch_index])
                nk.append(residualData[station_index]['ele'][epoch_index]+90)
        A = np.zeros((len(rk), 2))
        for i in range(0, len(rk)):
            A[i][0] = np.sin(np.deg2rad(nk[i])) * np.sin(np.deg2rad(ak[i]))
            A[i][1] = np.sin(np.deg2rad(nk[i])) * np.cos(np.deg2rad(ak[i]))

        AtA = A.transpose().dot(A)
        x = np.linalg.inv(AtA).dot(A.transpose()).dot(rk)
        yaw.append(np.rad2deg(np.arctan2(x[1], x[0])))
    return (epochs, yaw)


def solveLSEModel1WithRestraint(residualData):
    """

    :param residualData:
    :return:
    """
    epochs = []
    for station_index in residualData.keys():
        epochs = epochs + list(set(residualData[station_index]['mjd']) - set(epochs))
    epochs.sort()
    yaw = []

    for epoch in epochs:
        rk = []
        ak = []
        nk = []
        for station_index in residualData.keys():
            if epoch in residualData[station_index]['mjd']:
                epoch_index = residualData[station_index]['mjd'].index(epoch)

                rk.append(residualData[station_index]['rk'][epoch_index])
                ak.append(residualData[station_index]['azi'][epoch_index])
                nk.append(residualData[station_index]['ele'][epoch_index]+90)
        A = np.zeros((len(rk), 2))
        for i in range(0, len(rk)):
            A[i][0] = np.sin(np.deg2rad(nk[i])) * np.sin(np.deg2rad(ak[i]))
            A[i][1] = np.sin(np.deg2rad(nk[i])) * np.cos(np.deg2rad(ak[i]))


        AtA = A.transpose().dot(A)
        x = np.linalg.inv(AtA).dot(A.transpose()).dot(rk)
        y1 = np.rad2deg(np.arctan2(x[1], x[0]))
        y2 = y1
        for i in range(0, 50):
            A2 = np.append(A, [[2*constants.PCO_x*np.cos(np.deg2rad(y1)), 2*constants.PCO_x*np.sin(np.deg2rad(y1))]], axis=0)
            r2 = np.append(rk, [2*(constants.PCO_x**2)], axis=0)
            AtA = A2.transpose().dot(A2)
            x = np.linalg.inv(AtA).dot(A2.transpose()).dot(r2)
            y1 = y2
            y2 = np.rad2deg(np.arctan2(x[1], x[0]))
        yaw.append(y2)
    return (epochs, yaw)


def computeGrandDesignMatrixWindow(residualData, startEpoch, endEpoch):
    trimmedResData = utils.trimResDataWindow(residualData, startEpoch, endEpoch)
    (epochs_init, yaw_init) = solveLSEModel1WithRestraint(trimmedResData)

    Rk = np.empty((0, 0), float)
    grandDesignMatrixBlock1 = np.empty((0, 0), float)
    grandDesignMatrixBlock2 = np.empty((0, 0), float)
    grandDesignMatrixBlock3 = np.empty((0, 0), float)

    for epoch in epochs_init:
        rk = []
        ak = []
        nk = []
        for station_index in trimmedResData.keys():
            epoch_index = residualData[station_index]['mjd'].index(epoch)

            rk.append(residualData[station_index]['rk'][epoch_index])
            ak.append(residualData[station_index]['azi'][epoch_index])
            nk.append(residualData[station_index]['ele'][epoch_index]+90)
        A = np.zeros((len(rk), 2))
        for i in range(0, len(rk)):
            A[i][0] = np.sin(np.deg2rad(nk[i])) * np.sin(np.deg2rad(ak[i]))
            A[i][1] = np.sin(np.deg2rad(nk[i])) * np.cos(np.deg2rad(ak[i]))
        grandDesignMatrixBlock1 = block_diag(grandDesignMatrixBlock1, A)

        B = np.array([2*constants.PCO_x*np.cos(yaw_init[epochs_init.index(epoch)]),
                      2*constants.PCO_x*np.sin(yaw_init[epochs_init.index(epoch)])])
        grandDesignMatrixBlock2 = block_diag(grandDesignMatrixBlock2, B)

        C = np.zeros((len(rk), len(epochs_init)))
        C[:, epochs_init.index(epoch)] = 1
        C = np.concatenate((C, np.eye(len(rk))), axis=1)
        try:
            grandDesignMatrixBlock3 = np.concatenate((grandDesignMatrixBlock3, C), axis=0)
            Rk = np.concatenate((Rk, np.array(rk)), axis=0)
        except:
            grandDesignMatrixBlock3 = C
            Rk = np.array(rk)


    Rk = np.concatenate((Rk, 2*(constants.PCO_x**2)*np.ones(len(epochs_init))), axis=0)
    grandDesignMatrix = np.concatenate((grandDesignMatrixBlock1, grandDesignMatrixBlock2), axis=0)
    grandDesignMatrixBlock3 = np.concatenate((grandDesignMatrixBlock3, np.zeros((len(epochs_init), len(epochs_init)+len(rk)))), axis=0)
    grandDesignMatrix = np.concatenate((grandDesignMatrix, grandDesignMatrixBlock3), axis=1)

    return (grandDesignMatrix, Rk)

def solveLSEModel4WithRestraint(residualData, startEpoch, endEpoch):
    trimmedResData = utils.trimResDataWindow(residualData, startEpoch, endEpoch)
    (epochs, _) = solveLSEModel1WithRestraint(trimmedResData)

    (A, rk) = computeGrandDesignMatrixWindow(residualData, startEpoch, endEpoch)
    AtA = A.transpose().dot(A)
    X = np.linalg.inv(AtA).dot(A.transpose()).dot(rk)
    yaw = []
    for epoch in epochs:
        i = epochs.index(epoch)
        x = X[2*i: (2*i)+2]
        y = np.rad2deg(np.arctan2(x[1], x[0]))
        yaw.append(y)

    return (epochs, yaw)



def solveLSEModel4WithRestraintWindow(residualData, Ne=10):
    epochs = utils.getEpochsArray(residualData)

    epochs_final = []
    yaw = []
    for epoch in epochs:
        index = epochs.index(epoch)
        startEpoch = epoch
        try:
            endEpoch = epochs[index+Ne]
            e, y = solveLSEModel4WithRestraint(residualData, startEpoch=startEpoch, endEpoch=endEpoch)
            epochs_final.append(e[5])
            yaw.append(y[5])
        except:
            pass
    return (epochs_final, yaw)