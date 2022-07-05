import utils
import graphics
import estimation


if __name__ == '__main__':

    year = 21
    prn = 27
    doy = 59
    man = "N1"
    man_simulated = "{}_sim".format(man) # simulated .res files contain maneuvers with an S at the end to differentiate from real residuals

    # Get stations data, orbit data, stations list and residuals for considered maneuver, doy, and prn
    stationsData = utils.getStationsData()
    orbitData = utils.getOrbData(year, doy, prn, man)
    stationsList = utils.getStationsList(orbitData, stationsData)
    residualDataReal = utils.getResData(year, doy, prn, man)
    # residualDataReal = estimation.filterResiduals(residualDataReal)
    residualDataReal = utils.cleanResData(residualDataReal)


    # Create fake residual file, then read residuals and plot
    utils.simulatedResiduals(stationsList, orbitData, year, doy, prn, man_simulated)
    residualDataSimulated = utils.getResData(year, doy, prn, man_simulated)
    # graphics.plotResiduals(residualDataSimulated, stationsData, year, doy, prn, man_simulated)
    '''
    # Estimate yaw from simulated residuals, no restraint
    (epochs, yaw) = estimation.solveLSEModel1(residualDataSimulated, orbitData, stationsData)
    graphics.plotEstimatedYaw(epochs, yaw, orbitData, year, doy, prn, man_simulated)

    # Estimate yaw from real residuals, no restraint
    graphics.plotResiduals(residualDataReal, stationsData, year, doy, prn, man)
    (epochs, yaw) = estimation.solveLSEModel1(residualDataReal, orbitData, stationsData)
    graphics.plotEstimatedYaw(epochs, yaw, orbitData, year, doy, prn, man)
    '''


    # Estimate yaw from simulated residuals
    (epochs, yaw) = estimation.solveLSEModel1(residualDataSimulated, orbitData, stationsData)
    graphics.plotEstimatedYaw(epochs, yaw, orbitData, year, doy, prn, man, extension="_sim_model1")

    (epochs, yaw) = estimation.solveLSEModel2(residualDataSimulated, orbitData, stationsData)
    graphics.plotEstimatedYaw(epochs, yaw, orbitData, year, doy, prn, man, extension="_sim_model2")

    # Estimate yaw from real residuals
    graphics.plotResiduals(residualDataReal, stationsData, year, doy, prn, man)
    (epochs, yaw) = estimation.solveLSEModel1(residualDataReal, orbitData, stationsData)
    graphics.plotEstimatedYaw(epochs, yaw, orbitData, year, doy, prn, man, extension="_real_model1")

    (epochs, yaw) = estimation.solveLSEModel2(residualDataReal, orbitData, stationsData)
    graphics.plotEstimatedYaw(epochs, yaw, orbitData, year, doy, prn, man, extension="_real_model2")

    (epochs, yaw) = estimation.solveLSEModel3(residualDataReal, orbitData, stationsData)
    graphics.plotEstimatedYaw(epochs, yaw, orbitData, year, doy, prn, man, extension="_real_model3")

    (epochs, yaw) = estimation.solveLSEModel3(residualDataReal, orbitData, stationsData)
    graphics.plotEstimatedYaw(epochs, yaw, orbitData, year, doy, prn, man, extension="_real_model3")

    (epochs, yaw) = estimation.final4(residualDataReal, orbitData, stationsData, Ne=len(epochs)-1)
    graphics.plotEstimatedYaw(epochs, yaw, orbitData, year, doy, prn, man, extension="_real_model4")


