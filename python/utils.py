import os
import constants
import numpy as np
from estimation import getCLockCorrections


def getSlopeData(filename):
    slopeData = {}

    slopeData['name'] = []
    slopeData['inverted'] = []
    slopeData['delta_psi'] = []
    slopeData['t_delay'] = []
    slopeData['slope_n'] = []
    slopeData['slope_r'] = []
    slopeData['beta'] = []

    filepath = os.path.join(constants.OUT, filename)

    file = open(filepath)
    lines = file.readlines()
    for line in lines[1:]:
        items = [x for x in line.split('\t') if x != ""]
        slopeData['name'].append(items[0])

        slopeData['inverted'].append(int(items[1]))
        slopeData['delta_psi'].append(float(items[2]))
        slopeData['t_delay'].append(float(items[3]))

        slopeData['slope_n'].append(float(items[4]))
        slopeData['slope_r'].append(float(items[5]))
        slopeData['beta'].append(float(items[6]))

    return slopeData


def getCLockData(filename):
    clkData = {}

    clkData['mjd'] = []
    clkData['clk'] = []
    filepath = os.path.join(constants.OUT, filename)

    file = open(filepath)
    lines = file.readlines()
    for line in lines:
        items = [x for x in line.split(' ') if x != ""]
        clkData['mjd'].append(float(items[1]))
        clkData['clk'].append(float(items[2]))

    return clkData


def getOrbData(filename):
    """
    This function reads orbit files and stores the data in a dictionary names orbitData.
    This data is necessary to simulate measurements for yaw angles - and to compute nadir and azimuth angles.
    :param year: year as written in orb file name
    :param doy: day of year
    :param prn: prn of satellite (GPS block IIF)
    :param man: maneuver type - can be N1, N2, M1, M2
    :return: dictionary containing data read from .orb file
    """
    orbitData = {}

    orbitData['mjd'] = []
    orbitData['x'] = []
    orbitData['y'] = []
    orbitData['z'] = []
    orbitData['vx'] = []
    orbitData['vy'] = []
    orbitData['vz'] = []
    orbitData['beta'] = []
    orbitData['usun'] = []
    orbitData['yaw'] = []
    orbitData['shadow_flag'] = []

    filepath = os.path.join(constants.OUT, filename)

    file = open(filepath)
    lines = file.readlines()
    for line in lines:
        items = [x for x in line.split(' ') if x!=""]
        orbitData['mjd'].append(float(items[3]))

        orbitData['x'].append(float(items[4]))
        orbitData['y'].append(float(items[5]))
        orbitData['z'].append(float(items[6]))

        orbitData['vx'].append(float(items[7]))
        orbitData['vy'].append(float(items[8]))
        orbitData['vz'].append(float(items[9]))

        orbitData['beta'].append(float(items[10]))
        orbitData['usun'].append(float(items[11]))
        orbitData['yaw'].append(float(items[12]))
        orbitData['shadow_flag'].append(int(items[13]))

    return orbitData


def getSatellitePositionVelocity(orbitData, epoch):
    """
    This function takes an orbitData dictionary and an epoch and selects the position and velocity vectors of the
    satellite at that certain epoch
    :param orbitData: dictionary containing data read from .orb file
    :param epoch: epoch that we want to get position and velocity for
    :return: tuple containing position and velocity vectors as read from .orb file
    """
    for i in range(0, len(orbitData['mjd'])):
        if orbitData['mjd'][i] == epoch:
            r_sat = np.array([orbitData['x'][i], orbitData['y'][i], orbitData['z'][i]])
            v_sat = np.array([orbitData['vx'][i], orbitData['vy'][i], orbitData['vz'][i]])
    return (r_sat, v_sat)


def getResData(filename):
    """
    This function reads residual files and stores the data in a dictionary names residualData.
    This data is necessary to estimate yaw angles.
    :param year: year as written in orb file name
    :param doy: day of year
    :param prn: prn of satellite (GPS block IIF)
    :param man: maneuver type - can be N1, N2, M1, M2
    :return: dictionary containing data read from .res file
    """
    residualData = {}

    filepath = os.path.join(constants.OUT, filename)

    file = open(filepath)
    lines = file.readlines()
    for line in lines:
        items = [x for x in line.split(' ') if x!=""]
        station_index = int(items[0])
        if not station_index in residualData.keys():
            residualData[station_index] = {}
            residualData[station_index]['mjd'] = []
            residualData[station_index]['rk'] = []
            residualData[station_index]['azi'] = []
            residualData[station_index]['ele'] = []


        residualData[station_index]['mjd'].append(float(items[2]))
        residualData[station_index]['rk'].append(float(items[3]))
        residualData[station_index]['azi'].append(float(items[4]))
        residualData[station_index]['ele'].append(float(items[5]))

    return residualData



def correctResData(orbitData, stationsData, residualData, yaw_epochs, yaw):

    correctedData = {}
    stations = list(residualData.keys())

    for station_index in stations:
        (clock_epochs, clock_corrections, _) = getCLockCorrections(residualData, orbitData, stationsData)
        if not station_index in correctedData.keys():
            correctedData[station_index] = {}
            correctedData[station_index]['mjd'] = []
            correctedData[station_index]['rk'] = []
            correctedData[station_index]['mjd'] = residualData[station_index]['mjd']
            correctedData[station_index]['rk'] = [0 for i in range(0, len(residualData[station_index]['rk']))]

        for epoch in correctedData[station_index]['mjd']:
            sat_epoch_index = orbitData['mjd'].index(epoch)
            res_epoch_index = correctedData[station_index]['mjd'].index(epoch)
            yaw_epoch_index = yaw_epochs.index(epoch)
            clk_epoch_index = clock_epochs.index(epoch)

            r_sat = np.array([orbitData['x'][sat_epoch_index], orbitData['y'][sat_epoch_index], orbitData['z'][sat_epoch_index]])
            v_sat = np.array([orbitData['vx'][sat_epoch_index], orbitData['vy'][sat_epoch_index], orbitData['vz'][sat_epoch_index]])

            azi, nad = getAzimuthNadirSatellite(r_sat, v_sat, station_index, stationsData)
            rk = residual(azi, nad, yaw[yaw_epoch_index], noise=0)

            correctedData[station_index]['rk'][res_epoch_index] = residualData[station_index]['rk'][res_epoch_index] - rk - clock_corrections[clk_epoch_index]

    return correctedData



def cleanResData(residualData):
    """
    This function checks that the residual signal from a certain station does not have strange variations
    (as those that appear when many cycle slips are present). If a residual signal is detected which presents cycle slips
    or has a standard deviation above a certain threshold, then the signal from that station is removed.
    :param residualData: dictionary containing data from .res file
    :return: dictionary containing residual data with 'noisy' stations removed
    """
    stations = list(residualData.keys())
    for st in stations:
        residuals = list(residualData[st]['rk'])
        if np.std(residuals) > constants.CS_THRESHOLD:
            residualData.pop(st, None)
        for i in range(1, len(residuals)):
            if np.sqrt((residuals[i] - residuals[i-1])**2) > constants.CS_THRESHOLD:
                residualData.pop(st, None)

    return residualData


def trimResData(residualData, orbitData, stationsData):
    """
    This function selects only those stations that "see" a satellite for the whole duration of a maneuver.
    The residual data dictionary is trimmed in that the data from stations which only see the satellite for
    a part of the maneuver duration are removed.
    :param residualData: dictionary containing data from .res file
    :param orbitData: dictionary containing data from .orb file
    :param stationsData: dictionary containing stations info
    :return: trimmed residual data dictionary
    """
    stationsList = getStationsList(orbitData, stationsData)

    trimmedResData = {}
    for item in stationsList.keys():
        try:
            trimmedResData[item] = residualData[item]
        except KeyError:
            pass
    return trimmedResData


def getEpochsArray(residualData):
    """
    This function gets all unique epochs contained in a residual data dictionary.
    :param residualData: dictionary containing data from .res file
    :return: epochs list
    """
    epochs = []
    for station_index in residualData.keys():
        epochs = epochs + list(set(residualData[station_index]['mjd']) - set(epochs))
    epochs.sort()
    return epochs


def trimResDataWindow(residualData, startEpoch, endEpoch):
    """
    This function selects only those residual samples which correspond to the time interval between startEpoch and
    endEpoch. The data is rearranged into a new dictionary which is the return value.
    All selected stations contain data over all epochs so there are no gaps in data over any epochs.
    :param residualData: data read from .res file
    :param startEpoch: initial epoch of the considered time interval
    :param endEpoch: final epoch of the considered time interval
    :return: residual data dictionary trimmed between the two epochs
    """
    epochs = getEpochsArray(residualData)
    epochsList = epochs[epochs.index(startEpoch) : epochs.index(endEpoch)]
    trimmedResData = {}

    for station_index in residualData.keys():
        mjd = residualData[station_index]['mjd']
        if set(epochsList).issubset(mjd):
            indices = [mjd.index(item) for item in epochsList]

            trimmedResData[station_index] = {}
            trimmedResData[station_index]['mjd'] = []
            trimmedResData[station_index]['rk'] = []
            trimmedResData[station_index]['azi'] = []
            trimmedResData[station_index]['ele'] = []

            for index in indices:
                trimmedResData[station_index]['mjd'].append(residualData[station_index]['mjd'][index])
                trimmedResData[station_index]['rk'].append(residualData[station_index]['rk'][index])
                trimmedResData[station_index]['azi'].append(residualData[station_index]['azi'][index])
                trimmedResData[station_index]['ele'].append(residualData[station_index]['ele'][index])
    return trimmedResData


def getLonLat(x, y, z):
    """
    This function gets a position in x, y, z (inertial frame - ECEF) and returns the longitude and latitude (in radians) in a tuple.
    """
    lon = np.arctan2(y, x)
    lat = np.arctan2(z, np.sqrt(x**2 + y**2))
    return (lon, lat)


def getAzimuthNadirSatellite(r_sat, v_sat, stationIndex, stationsData):
    """
    This function computed the azimuth and nadir angles under which a satellite is observed by a certain station.
    The angles are in satellite fixed frame.
    :param r_sat: satellite position vector
    :param v_sat: satellite velocity vector
    :param stationIndex: index of station considered
    :param stationsData: dictionary containing data for all stations (such as index, name and position)
    :return:
    """
    x_station = stationsData[stationIndex]['x']
    y_station = stationsData[stationIndex]['y']
    z_station = stationsData[stationIndex]['z']

    r_st = np.array([x_station, y_station, z_station])
    r = np.sqrt(r_sat[0]**2 + r_sat[1]**2 + r_sat[2]**2)

    h = np.cross(r_sat, v_sat)
    gamma = -(1/r**2)*(r_st.dot(r_sat))
    R_proj = r_st + gamma*r_sat
    R_r = r_st - r_sat
    nad_nominator = np.sqrt(R_r[0]**2 + R_r[1]**2 + R_r[2]**2)*r

    azi = np.arctan2(R_proj.dot(v_sat), (1/r)*(R_proj.dot(h)))
    nad = np.arccos(-(R_r.dot(r_sat))/nad_nominator)

    return (azi, nad)


def getAzimuthElevationTopocentric(r_sat, stationIndex, stationsData):
    """
    This function computes the azimuth and elevation of a satellite from a given station.
    :param x_sat: satellite x position - ECEF
    :param y_sat: satellite y position - ECEF
    :param z_sat: satellite z position - ECEF
    :param stationIndex: index of wanted station - ECEF
    :param stationsData: dictionary given by getStationsList - contains all stations positions, name and indices
    :return: azimuth and elevation angles in radians, grouped in a tuple
    """
    x_station = stationsData[stationIndex]['x']
    y_station = stationsData[stationIndex]['y']
    z_station = stationsData[stationIndex]['z']

    r_station = np.array([x_station, y_station, z_station])
    r_trans = r_sat - r_station
    lon, lat = getLonLat(x_station, y_station, z_station)

    Q1 = np.array([[-1, 0, 0], [0, 1, 0], [0, 0, 1]])
    R2 = np.array([[np.cos((np.pi/2)-lat), 0, -np.sin((np.pi/2)-lat)], [0, 1, 0], [np.sin((np.pi/2)-lat), 0, np.cos((np.pi/2)-lat)]]) # R2((np.pi/2)-lat)
    R3 = np.array([[np.cos(lon), np.sin(lon), 0], [-np.sin(lon), np.cos(lon), 0], [0, 0, 1]]) # R3(lon)
    r4 = Q1.dot(R2).dot(R3).dot(r_trans)

    a = np.arctan2(r4[1], r4[0]) + np.pi*2
    if np.rad2deg(a) >= 360:
        a = np.arctan2(r4[1], r4[0])

    e = np.arctan2(r4[2], np.sqrt(r4[0]**2 + r4[1]**2))
    return (a, e)


def isVisible(r_sat, stationIndex, stationsData):
    """
    This function tells us whether a satellite is visible from a certain station based on satellite position at a certain epoch.
    :param x_sat: satellite x position - ECEF
    :param y_sat: satellite y position - ECEF
    :param z_sat: satellite z position - ECEF
    :param stationIndex: index of wanted station - ECEF
    :param stationsData: dictionary given by getStationsList - contains all stations positions, name and indices
    :return: True if satellite elevation is above constants.ELEVATION_MASK and False otherwise
    """
    if getAzimuthElevationTopocentric(r_sat, stationIndex, stationsData)[1] >= np.deg2rad(constants.ELEVATION_MASK):
        return True
    else:
        return False


def getStationsData():
    """
    This function reads all stations from the file GPSYAW.crd.
    :return: stationsData is a dictionary containing station indices as keys, and their positions (x, y, z and the name) as values
    """
    stationsData = {}

    file = open(constants.STATION_COORDINATE_LIST)
    lines = file.readlines()
    for line in lines:
        items = [x for x in line.split(' ') if x != ""]
        if len(items) > 3:
            stationsData[int(items[0])] = {}
            stationsData[int(items[0])]['x'] = float(items[1])
            stationsData[int(items[0])]['y'] = float(items[2])
            stationsData[int(items[0])]['z'] = float(items[3])
            stationsData[int(items[0])]['name'] = items[4].strip()
    return stationsData


def getStationsList(orbitData, stationsData):
    """
    This function computes an intersection between .orb and .crd files.
    We want to see which stations see the satellite for the whole duration of the maneuver.
    :param orbitData: dictionary containing data read from .orb file, created by getOrbData()
    :param stationsData: dictionary with stations coordinates created by getStationsData()
    :return: dictionary containing list of stations with index and name, fitted to given .orb file
    """
    stations_list = {}
    for station_index in stationsData.keys():
        r_sat_init = np.array([orbitData['x'][0], orbitData['y'][0],  orbitData['z'][0]])
        r_sat_end = np.array([orbitData['x'][-1], orbitData['y'][-1], orbitData['z'][-1]])
        if (isVisible(r_sat_init, stationIndex=station_index, stationsData=stationsData) and
            isVisible(r_sat_end, stationIndex=station_index, stationsData=stationsData)):
            stations_list[station_index] = {}
            stations_list[station_index]['name'] = stationsData[station_index]['name']
            stations_list[station_index]['x'] = stationsData[station_index]['x']
            stations_list[station_index]['y'] = stationsData[station_index]['y']
            stations_list[station_index]['z'] = stationsData[station_index]['z']
    return stations_list


def residual(azi, nad, yaw, noise=0.001):
    """
    This function computes a residual sample at a given epoch for a given station.
    :param azi: azimuth angle in radians
    :param nad: nadir angle in radians
    :param yaw: nominal yaw angle in degrees
    :param noise: maximum value of added noise in meters
    :return: residual sample rk in meters
    """
    yaw = np.deg2rad(yaw)
    rk = constants.PCO_x * np.cos(yaw) * np.sin(azi) * np.sin(nad) + \
         constants.PCO_x * np.sin(yaw) * np.cos(azi) * np.sin(nad) + \
         np.random.normal(0, noise, 1)[0]
    return rk


def simulatedResiduals(stationsList, orbitData, year=21, doy=58, prn=27, man="M1S", noise=0.001):
    """
    This function reates a fake residual log file, in the exact same format as the .res files created by BPE.
    This format is used in order to minimize the written code - the same functions which are used for
    real .res files work on simulated residuals as well.
    :param stationsList: dictionary containing information on stations which see the satellite for the whole duration of the maneuver
    :param orbitData: dictionary containing data read from .orb file
    :param prn: satellite PRN code
    :param filename: name given to fake .res file
    :return: writes .res file, with no return value
    """
    filename = "YAW" + str(year) + "0" + str(doy) + "0" + man + "G" + str(prn) + "_sim.res"
    filepath = os.path.join(constants.OUT, filename)
    file = open(filepath, "w")

    for i in range(0, len(orbitData['x'])):
        for station_index in stationsList.keys():
            r_sat = np.array([orbitData['x'][i], orbitData['y'][i], orbitData['z'][i]])
            v_sat = np.array([orbitData['vx'][i], orbitData['vy'][i], orbitData['vz'][i]])

            azi, nad = getAzimuthNadirSatellite(r_sat, v_sat, station_index, stationsList)
            _, ele = getAzimuthElevationTopocentric(r_sat, station_index, stationsList)

            rk = residual(azi, nad, orbitData['yaw'][i], noise=noise)
            mjd = orbitData['mjd'][i]
            line = "{} {} {} {} {} {}\n".format(station_index, prn, mjd, rk, np.rad2deg(azi), np.rad2deg(ele))
            file.write(line)
    file.close()


def logInfo(yaw, yaw_n, errors, filename, extension):
    yaw = np.array(yaw)
    yaw_n = np.array(yaw_n)
    errors = np.array(errors)

    std_dev = np.std(yaw-yaw_n)
    rms = np.sqrt(np.mean(np.square(yaw-yaw_n)))
    err_mean = np.mean(errors)

    logName = filename[0:14] + "_{}.txt".format(extension)
    subdir = filename[0:14]
    if not os.path.exists(os.path.join(constants.FIGS, subdir)):
        os.makedirs(os.path.join(constants.FIGS, subdir))
    logPath = os.path.join(constants.FIGS, subdir, logName)

    file = open(logPath, "a")
    file.write("{}\t {}\t {}\n".format(std_dev, rms, err_mean))
    file.close()





