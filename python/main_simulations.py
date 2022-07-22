import utils
import graphics
import estimation
import constants
import os

if __name__ == '__main__':

    stationsData = utils.getStationsData()

    file = "YAW210590M1G27.res"
    year = 21
    doy = 59
    prn = 27
    man = file[9:11]
    man_simulated = "{}_sim".format(man)

    # Get stations data, orbit data, stations list and residuals for considered maneuver, doy, and prn
    orbfile = str(file[0:14]) + ".orb"
    orbitData = utils.getOrbData(orbfile)
    stationsList = utils.getStationsList(orbitData, stationsData)
    Ne = len(orbitData['mjd']) - 1

    noiseValues = [x/10000 for x in range(1, 200, 5)]
    for noise in noiseValues:

        # Create fake residual file, then read residuals and plot
        utils.simulatedResiduals(stationsList, orbitData, year, doy, prn, man, noise=noise)
        residualDataSimulated = utils.getResData("YAW210590M1G27_sim.res")
        graphics.plotResiduals(residualDataSimulated, stationsData, file)

        # Estimate yaw from simulated residuals, no restraint
        (epochs, yaw, errors) = estimation.solveLSEModel1(residualDataSimulated, orbitData, stationsData)
        graphics.plotEstimatedYawErrorbars(epochs, yaw, errors, orbitData, file,
                                           extension="_sim_model1_bars")
        graphics.plotEstimatedYaw(epochs, yaw, orbitData, file,
                                           extension="_sim_model1")
        utils.logInfo(yaw, orbitData['yaw'], errors, file, extension="info_model1")


        (epochs, yaw, errors) = estimation.solveLSEModel2(residualDataSimulated, orbitData, stationsData)
        graphics.plotEstimatedYawErrorbars(epochs, yaw, errors, orbitData, file,
                                           extension="_sim_model2_bars")
        graphics.plotEstimatedYaw(epochs, yaw, orbitData, file,
                                  extension="_sim_model2")
        utils.logInfo(yaw, orbitData['yaw'], errors, file, extension="info_model2")


        (epochs, yaw, errors) = estimation.solveLSEModel3(residualDataSimulated, orbitData, stationsData)
        graphics.plotEstimatedYawErrorbars(epochs, yaw, errors, orbitData, file,
                                           extension="_sim_model3_bars")
        graphics.plotEstimatedYaw(epochs, yaw, orbitData, file,
                                  extension="_sim_model3")
        utils.logInfo(yaw, orbitData['yaw'], errors, file, extension="info_model3")


        (epochs, yaw, errors) = estimation.solveLSEModel4(residualDataSimulated, orbitData, stationsData, Ne=Ne)
        graphics.plotEstimatedYawErrorbars(epochs, yaw, errors, orbitData, file,
                                           extension="_sim_model4_bars")
        graphics.plotEstimatedYaw(epochs, yaw, orbitData, file,
                                  extension="_sim_model4")
        utils.logInfo(yaw, orbitData['yaw'][0:len(yaw)], errors, file, extension="info_model4")

    graphics.plotLogInfo(noiseValues, file, extension="info_model1")
    graphics.plotLogInfo(noiseValues, file, extension="info_model2")
    graphics.plotLogInfo(noiseValues, file, extension="info_model3")
    graphics.plotLogInfo(noiseValues, file, extension="info_model4")

