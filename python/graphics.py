import constants
import matplotlib.pyplot as plt
import os


def plotResiduals(residualData, stationsData, year=21, doy=58, prn=27, man="M1"):
    """
    This function plots the data read from .res files and saves them in a .png figure.
    :param residualData: dictionary containing residual data as created by utils.getResData()
    :param stationsData: dictionary containing station data as created by utils.getStationsData()
    :param year: year - used to construct .res file name
    :param prn: satellite prn code - used to construct .res file name
    :param doy: day of year - used to construct .res file name
    :param man: maneuver type - used to construct .res file name
    :return: writes .jpg figure, with no return value
    """
    figName = str(year) + "0" + str(doy) + "0" + man + "G" + str(prn) + ".jpg"
    figPath = os.path.join(constants.FIGS, figName)
    plt.figure(figsize=(12, 6))

    for station_index in residualData.keys():
        tk = residualData[station_index]['mjd']
        rk = residualData[station_index]['rk']
        st_name = stationsData[station_index]['name']
        plt.plot(tk, rk, label=st_name)
    plt.legend(loc="upper left")
    plt.grid()
    plt.xlabel("Time (MJD)")
    plt.ylabel("Residual (meters)")
    plt.title("Residual values for PRN " + str(prn) + ", DOY " + str(doy) + ", maneuver type " + man)
    plt.savefig(figPath)


def plotNominalYaw(orbitData, year=21, prn=27, doy=58, man="M1", savefig=True):
    """
    This function plots the nominal yaw angle read from .orb files, alongside the shadow marking.
    :param orbitData: dictionary containing data read from .orb file, as created by utils.getOrbData()
    :param year: year - used to construct .res figure name
    :param prn: satellite prn code - used to construct figure name
    :param doy: day of year - used to construct figure name
    :param man: maneuver type - used to construct figure name
    :return: writes .jpg figure, with no return value
    """
    figName = str(year) + "0" + str(doy) + "0" + man + "G" + str(prn) + "_yawn.jpg"
    figPath = os.path.join(constants.FIGS, figName)
    plt.figure(figsize=(12, 6))

    time = orbitData['mjd']
    yaw = orbitData['yaw']
    shadow = orbitData['shadow_flag']
    beta = [str(orbitData['beta'][0]), str(orbitData['beta'][-1])]

    plt.plot(time, yaw)
    for i in range(0, len(shadow)-1):
        if shadow[i] == 1:
            plt.axvspan(time[i], time[i+1], color='lightgrey', alpha=0.75, lw=0, zorder=-1)
    plt.grid()
    plt.xlabel("Time (MJD)")
    plt.ylabel("Yaw angle (degrees)")
    plt.title("Yaw angle values for PRN " + str(prn) + ", DOY " + str(doy) + ", maneuver type " +
              man + ", ß: " + beta[0] + "° - " + beta[1] + "°")
    if savefig:
        plt.savefig(figPath)
    else:
        return plt


def plotEstimatedYaw(epochs, yaw, orbitData, year=21, doy=58, prn=27, man="M1"):
    figName = str(year) + "0" + str(doy) + "0" + man + "G" + str(prn) + "_yawest.jpg"
    figPath = os.path.join(constants.FIGS, figName)
    plt = plotNominalYaw(orbitData, year, prn, doy, man, savefig=False)
    plt.scatter(epochs, yaw, marker='.', color='r', zorder=1)
    plt.savefig(figPath)