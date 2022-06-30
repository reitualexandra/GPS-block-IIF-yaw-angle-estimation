import utils
import graphics
import estimation


if __name__ == '__main__':

    year = 21
    prn = 27
    doy = 58
    man = "M1"
    man_simulated = "M1S" # simulated .res files contain maneuvers with an S at the end to differentiate from real residuals
    man_restraint = "M1R"
    man_simulated_restraint = "M1RS"

    # Get stations data, orbit data, stations list and residuals for considered maneuver, doy, and prn
    stationsData = utils.getStationsData()
    orbitData = utils.getOrbData(year, doy, prn, man)
    stationsList = utils.getStationsList(orbitData, stationsData)
    residualDataReal = utils.getResData(year, doy, prn, man)

    # Create fake residual file, then read residuals and plot
    utils.simulatedResiduals(stationsList, orbitData, year, doy, prn, man_simulated)
    residualDataSimulated = utils.getResData(year, doy, prn, man_simulated)
    graphics.plotResiduals(residualDataSimulated, stationsData, year, doy, prn, man_simulated)

    # Estimate yaw from simulated residuals, no restraint
    (epochs, yaw) = estimation.solveLSEModel1(residualDataSimulated)
    graphics.plotEstimatedYaw(epochs, yaw, orbitData, year, doy, prn, man_simulated)

    # Estimate yaw from real residuals, no restraint
    graphics.plotResiduals(residualDataReal, stationsData, year, doy, prn, man)
    (epochs, yaw) = estimation.solveLSEModel1(residualDataReal)
    graphics.plotEstimatedYaw(epochs, yaw, orbitData, year, doy, prn, man)

    # Estimate yaw from real residuals, with restraint
    # residualDataReal = utils.trimResData(residualDataReal, orbitData, stationsData)
    (epochs, yaw) = estimation.solveLSEModel1WithRestraint(residualDataReal)
    graphics.plotEstimatedYaw(epochs, yaw, orbitData, year, doy, prn, man_restraint)

    # Estimate yaw from simulated residuals, with restraint
    (epochs, yaw) = estimation.solveLSEModel1WithRestraint(residualDataSimulated)
    graphics.plotEstimatedYaw(epochs, yaw, orbitData, year, doy, prn, man_simulated_restraint)
