import constants
import estimation
import utils
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 16})
import numpy as np
import os
import datetime


def mjd2hms(mjd):
    times = []

    for item in mjd:
        hms = 10**5 * (item - int(item))
        hh = str(datetime.timedelta(seconds=hms))[0:5]
        if hh[-1] == ":":
            hh = hh[0:4]
        times.append(hh)
    return times


def getDetails(filename):
    year = filename[3:5]
    doy = filename[5:9]
    prn = filename[11:14]
    man = filename[9:11]
    return year, doy, man, prn


def plotResiduals(residualData, stationsData, filename):
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
    year, doy, man, prn = getDetails(filename)
    figName = filename[0:14] + ".jpg"

    subdir = filename[0:14]
    if not os.path.exists(os.path.join(constants.FIGS, subdir)):
        os.makedirs(os.path.join(constants.FIGS, subdir))

    figPath = os.path.join(constants.FIGS, subdir, figName)
    plt.figure(figsize=(12, 6))

    for station_index in residualData.keys():
        tk = residualData[station_index]['mjd']
        rk = residualData[station_index]['rk']
        st_name = stationsData[station_index]['name']

        plt.plot(tk, rk, label=st_name)

    tk = utils.getEpochsArray(residualData)
    labels = mjd2hms(tk)
    plt.xticks(tk[::60], labels[::60])

    #plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), mode="expand")
    plt.grid()
    plt.xlabel("Time")
    plt.ylabel("Residual (m)")
    plt.title("PRN " + str(prn) + ", DOY " + str(doy)[0:3] + ", " + man)
    plt.savefig(figPath)
    #plt.show()
    plt.close('all')


def plotResidualsCorrected(residualData, correctedResidualData, filename, extension=""):
    year, doy, man, prn = getDetails(filename)
    figName = filename[0:14] + "{}.jpg".format(extension)

    subdir = filename[0:14]
    if not os.path.exists(os.path.join(constants.FIGS, subdir)):
        os.makedirs(os.path.join(constants.FIGS, subdir))

    figPath = os.path.join(constants.FIGS, subdir, figName)
    plt.figure(figsize=(12, 6))

    for station_index in residualData.keys():
        tk1 = residualData[station_index]['mjd']
        rk1 = residualData[station_index]['rk']
        plt.scatter(tk1, rk1, marker='^', color='lime')

    for station_index in residualData.keys():
        tk2 = correctedResidualData[station_index]['mjd']
        rk2 = correctedResidualData[station_index]['rk']
        plt.scatter(tk2, rk2, marker='.', color='fuchsia')

    tk = utils.getEpochsArray(residualData)
    labels = mjd2hms(tk)
    plt.xticks(tk[::60], labels[::60])

    lime_patch = mpatches.Patch(color='lime', label='Original residuals')
    purple_patch = mpatches.Patch(color='fuchsia', label='Corrected residuals')
    plt.legend(handles=[lime_patch, purple_patch], loc='lower right')
    # plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), mode="expand")
    plt.grid()
    plt.xlabel("Time")
    plt.ylabel("Residual (m)")
    plt.title("PRN " + str(prn) + ", DOY " + str(doy)[0:3] + ", " + man)
    plt.savefig(figPath)
    #plt.show()
    plt.close('all')



def plotNominalYaw(orbitData, filename, savefig=False):
    """
    This function plots the nominal yaw angle read from .orb files, alongside the shadow marking.
    :param orbitData: dictionary containing data read from .orb file, as created by utils.getOrbData()
    :param year: year - used to construct .res figure name
    :param doy: day of year - used to construct figure name
    :param prn: satellite prn code - used to construct figure name
    :param man: maneuver type - used to construct figure name
    :return: writes .jpg figure, with no return value
    """
    year, doy, man, prn = getDetails(filename)
    figName = filename[0:14] + "_yawn.jpg"

    subdir = filename[0:14]
    if not os.path.exists(os.path.join(constants.FIGS, subdir)):
        os.makedirs(os.path.join(constants.FIGS, subdir))

    figPath = os.path.join(constants.FIGS, subdir, figName)

    time = orbitData['mjd']
    yaw = orbitData['yaw']
    shadow = orbitData['shadow_flag']
    beta = [str(orbitData['beta'][0]), str(orbitData['beta'][-1])]

    plt.figure(figsize=(12, 6))
    plt.plot(time, yaw, linestyle='dashed', color="red", label="Modeled")

    labels = mjd2hms(time)
    plt.xticks(time[::60], labels[::60])
    #plt.xticks(time, labels)

    for i in range(0, len(shadow)-1):
        if shadow[i] == 1:
            plt.axvspan(time[i], time[i+1], color='lightgrey', alpha=0.75, lw=0, zorder=-1)
    plt.grid()
    plt.xlabel("Time")
    plt.ylabel("Ψ angle (deg)")
    plt.title("PRN " + str(prn) + ", DOY " + str(doy)[0:3] + ", " +
              man + ", " + min(beta) + "° < ß < " + max(beta) + "°")
    if savefig:
        plt.savefig(figPath)
        plt.close('all')
    else:
        return plt


def plotEstimatedYaw(epochs, yaw, orbitData, filename, extension=""):
    """
    This function plots the estimated yaw values over the nominal ones read from .orb data file.
    :param epochs: epochs over which yaw was estimated
    :param yaw: estimated yaw values in degrees
    :param orbitData: dictionary containing data read from .orb file, as created by utils.getOrbData()
    :param year: year - used to construct .res figure name
    :param doy: day of year - used to construct figure name
    :param prn: satellite prn code - used to construct figure name
    :param man: maneuver type - used to construct figure name
    :return:
    """
    figName = filename[0:14] + "_yawest{}.jpg".format(extension)
    yaw = [x * (-1) for x in yaw]

    subdir = filename[0:14]
    if not os.path.exists(os.path.join(constants.FIGS, subdir)):
        os.makedirs(os.path.join(constants.FIGS, subdir))

    figPath = os.path.join(constants.FIGS, subdir, figName)
    plt = plotNominalYaw(orbitData, filename)
    plt.scatter(epochs, yaw, marker='.', color='cornflowerblue', zorder=1, label="Measured")

    plt.legend(loc="upper right")
    plt.savefig(figPath)
    plt.close('all')



def plotEstimatedYawInterpolated(epochs, yaw, orbitData, filename, extension=""):
    """
    This function plots the estimated yaw values over the nominal ones read from .orb data file.
    :param epochs: epochs over which yaw was estimated
    :param yaw: estimated yaw values in degrees
    :param orbitData: dictionary containing data read from .orb file, as created by utils.getOrbData()
    :param year: year - used to construct .res figure name
    :param doy: day of year - used to construct figure name
    :param prn: satellite prn code - used to construct figure name
    :param man: maneuver type - used to construct figure name
    :return:
    """
    figName = filename[0:14] + "_yawest{}.jpg".format(extension)
    yaw = [x * (-1) for x in yaw]

    subdir = filename[0:14]
    if not os.path.exists(os.path.join(constants.FIGS, subdir)):
        os.makedirs(os.path.join(constants.FIGS, subdir))

    figPath = os.path.join(constants.FIGS, subdir, figName)
    plt = plotNominalYaw(orbitData, filename)
    plt.scatter(epochs, yaw, marker='.', color='cornflowerblue', zorder=1, label="Measured")

    y_filt = estimation.interpolate(epochs, yaw)
    plt.plot(epochs, y_filt, color='tomato')

    plt.legend(loc="upper right")
    plt.savefig(figPath)
    plt.close('all')


def plotEstimatedYawErrorbars(epochs, yaw, errors, orbitData, filename, extension=""):
    """
    This function plots the estimated yaw values over the nominal ones read from .orb data file.
    :param epochs: epochs over which yaw was estimated
    :param yaw: estimated yaw values in degrees
    :param orbitData: dictionary containing data read from .orb file, as created by utils.getOrbData()
    :param year: year - used to construct .res figure name
    :param doy: day of year - used to construct figure name
    :param prn: satellite prn code - used to construct figure name
    :param man: maneuver type - used to construct figure name
    :return:
    """
    figName = filename[0:14] + "_yawest{}.jpg".format(extension)
    #errors = [x*0.5 for x in errors]
    yaw = [x * (-1) for x in yaw]

    subdir = filename[0:14]
    if not os.path.exists(os.path.join(constants.FIGS, subdir)):
        os.makedirs(os.path.join(constants.FIGS, subdir))

    figPath = os.path.join(constants.FIGS, subdir, figName)
    plt = plotNominalYaw(orbitData, filename, savefig=False)

    plt.errorbar(epochs, yaw, yerr=errors, ecolor="powderblue", zorder=0, ls='none')
    plt.scatter(epochs, yaw, marker='.', color='cornflowerblue', zorder=1, label="Measured")

    plt.legend(loc="upper right")
    plt.savefig(figPath)
    #plt.show()
    plt.close('all')



def plotClockCorrections(clk, clkData, filename, extension=""):
    figName = filename[0:14] + "_clk{}.jpg".format(extension)

    subdir = filename[0:14]
    if not os.path.exists(os.path.join(constants.FIGS, subdir)):
        os.makedirs(os.path.join(constants.FIGS, subdir))

    figPath = os.path.join(constants.FIGS, subdir, figName)

    time = clkData['mjd']
    clk_n = clkData['clk']

    t = np.array(time)
    m, b = np.polyfit(t, clk_n, 1)
    clk_interpolated = m * t + b
    clk_n = clk_n - clk_interpolated

    plt.figure(figsize=(12, 6))
    plt.plot(time, clk_n, marker='.', color="red", label="Old clock correction")

    clk = np.divide(clk, constants.speed_of_light)
    plt.plot(time, clk_n + clk, marker='.', color="cornflowerblue", label="New clock correction")
    #plt.plot(time, clk_n - clk, marker='.', color="yellow", label="SECOND New clock correction")

    labels = mjd2hms(time)
    plt.xticks(time[::60], labels[::60])

    plt.grid()
    plt.xlabel("Time")
    plt.ylabel("Clock correction (sec)")

    plt.legend(loc="upper right")
    plt.savefig(figPath)
    #plt.show()
    plt.close('all')



def plotLogInfo(noiseValues, filename, extension):
    logName = filename[0:14] + "_{}.txt".format(extension)
    figName = filename[0:14] + "_{}".format(extension)
    subdir = filename[0:14]
    if not os.path.exists(os.path.join(constants.FIGS, subdir)):
        os.makedirs(os.path.join(constants.FIGS, subdir))
    logPath = os.path.join(constants.FIGS, subdir, logName)
    figPath = os.path.join(constants.FIGS, subdir, figName)

    file = open(logPath, "r")
    lines = file.readlines()

    stds = []
    dists= []
    errs = []

    for line in lines:
        items = [x for x in line.split(' ') if x != ""]
        stds.append(float(items[0]))
        dists.append(float(items[1]))
        errs.append(float(items[2]))

    noiseValues = np.array(noiseValues)*1000

    plt.figure(figsize=(12, 6))
    plt.scatter(noiseValues, stds, marker='o', color='cornflowerblue', zorder=1)

    m, b = np.polyfit(noiseValues, stds, 1)
    plt.plot(noiseValues, m * noiseValues + b, linestyle='dashed', color="red")

    plt.grid()
    plt.xlabel("Noise (mm)")
    plt.ylabel("STD (deg)")
    plt.title("Standard deviation variation with noise")
    plt.savefig(figPath + "_std.jpg")

    plt.figure(figsize=(12, 6))
    plt.scatter(noiseValues, dists, marker='o', color='cornflowerblue', zorder=1)

    m, b = np.polyfit(noiseValues, dists, 1)
    plt.plot(noiseValues, m * noiseValues + b, linestyle='dashed', color="red")

    plt.grid()
    plt.xlabel("Noise (mm)")
    plt.ylabel("RMS (deg)")
    plt.title("Root mean square variation with noise")
    plt.savefig(figPath + "_rms.jpg")

    plt.figure(figsize=(12, 6))
    plt.scatter(noiseValues, errs, marker='o', color='cornflowerblue', zorder=1)

    m, b = np.polyfit(noiseValues, errs, 1)
    plt.plot(noiseValues, m * noiseValues + b, linestyle='dashed', color="red")

    plt.grid()
    plt.xlabel("Noise (mm)")
    plt.ylabel("Mean errors (deg)")
    plt.title("Mean formal errors variation with noise")
    plt.savefig(figPath + "_err.jpg")
    plt.close('all')


def plotSlopeOverBeta(slopeData, extension=""):
    figName = "SlopeOverBeta_{}.jpg".format(extension)
    figPath = os.path.join(constants.FIGS, figName)

    plt.figure(figsize=(12, 6))
    plt.scatter(slopeData['beta'], slopeData['slope_n'], marker='o', color='fuchsia', zorder=1, label="Modeled")
    plt.scatter(slopeData['beta'], slopeData['slope_r'], marker='^', color='lime', zorder=1, label="Measured")

    plt.grid()
    plt.xlabel("ß (deg)")
    plt.ylabel("Slope (deg/s)")
    plt.title("Slope of nominal model over ß angle\n Estimated median slope m={:.3f} (deg/s)".format(np.median(slopeData['slope_r'])))
    plt.legend(loc="upper right")
    plt.savefig(figPath)



def plotInvertedOverBeta(slopeData, extension=""):
    figName = "SlopeInvOverBeta_{}.jpg".format(extension)
    figPath = os.path.join(constants.FIGS, figName)

    plt.figure(figsize=(8, 6))

    invBetas = []
    invVal = []
    normBetas = []
    normVal = []
    for i in range(0, len(slopeData['beta'])):
        if slopeData['inverted'][i]==1:
            invBetas.append(slopeData['beta'][i])
            invVal.append(slopeData['inverted'][i])
        else:
            normBetas.append(slopeData['beta'][i])
            normVal.append(slopeData['inverted'][i])

    plt.scatter(invBetas, invVal, marker='o', color='fuchsia', zorder=1, label="Inverted", s=80)
    plt.scatter(normBetas, normVal, marker='^', color='lime', zorder=1, label="Normal", s=80)
    plt.yticks([-0.25, 0, 0.5, 1, 1.25], ["", "+", "", "-", ""])

    plt.grid()
    plt.xlabel("ß (deg)")
    plt.ylabel("Slope sign")
    plt.title("Slope inversion over ß angle")
    plt.legend(loc="upper right")
    plt.savefig(figPath)


def plotEstimatedSlopeOverBeta(slopeData, extension=""):
    figName = "SlopeEstOverBeta_{}.jpg".format(extension)
    figPath = os.path.join(constants.FIGS, figName)

    noonSlopes = []
    noonBetas = []
    midnSlopes = []
    midnBetas = []
    for name in slopeData['name']:
        index = slopeData['name'].index(name)
        if name[9]=="N":
            noonSlopes.append(slopeData['slope_r'][index])
            noonBetas.append(slopeData['beta'][index])
        else:
            midnSlopes.append(slopeData['slope_r'][index])
            midnBetas.append(slopeData['beta'][index])


    plt.figure(figsize=(12, 6))
    plt.scatter(noonBetas, noonSlopes, marker='o', color='fuchsia', zorder=1, label="Noon turns")
    plt.scatter(midnBetas, midnSlopes, marker='^', color='lime', zorder=1, label="Midnight turns")

    plt.grid()
    plt.xlabel("ß (deg)")
    plt.ylabel("Slope (deg/s)")
    plt.title("Estimated slope over ß angle. Estimated median slope:\n noon: {:.3f} (deg/s), midn: {:.3f} (deg/s)".format(
        np.median(noonSlopes), np.median(midnSlopes)))
    plt.legend(loc="upper right")
    plt.savefig(figPath)


