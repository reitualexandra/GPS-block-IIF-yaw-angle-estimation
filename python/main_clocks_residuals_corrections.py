import utils
import graphics
import estimation
import constants
import os

if __name__ == '__main__':

    stationsData = utils.getStationsData()

    files = os.listdir(constants.OUT)
    res_files = [x for x in files if (".res" in x or ".RES" in x) and x.startswith("YAW")]
    for file in res_files: #YAW203580M1G25.RES

        try:
            file = str(file)
            year = file[3:5]
            doy = file[5:9]
            prn = file[11:14]
            man = file[9:11]
            man_simulated = "{}_sim".format(man)

            # Get stations data, orbit data, stations list and residuals for considered maneuver, doy, and prn
            orbfile = "YAW" + str(file[3:14]) + ".orb"  #YAW203580M1G25.orb
            clkFile = str(file[0:14]) + ".clk"

            orbitData = utils.getOrbData(orbfile)
            clockData = utils.getCLockData(clkFile)
            Ne = len(orbitData['mjd']) - 1

            stationsList = utils.getStationsList(orbitData, stationsData)
            residualDataReal = utils.getResData(file)
            residualDataReal = estimation.filterResiduals(residualDataReal)
            residualDataReal = utils.cleanResData(residualDataReal)

            graphics.plotResiduals(residualDataReal, stationsData, file)

            #(epochs, clk, bias) = estimation.solveLSEModel4ClocksBiases(residualDataReal, orbitData, stationsData, Ne=Ne)

            (epochs, clock_corrections, errors) = estimation.getCLockCorrections(residualDataReal, orbitData, stationsData)
            graphics.plotClockCorrections(clock_corrections, clockData, file, extension="_clocks")
            '''
            (epochs, yaw, errors) = estimation.solveLSEModel1(residualDataReal, orbitData, stationsData)
            correctedResiduals = utils.correctResData(orbitData, stationsData, residualDataReal, epochs, yaw)
            graphics. plotResidualsCorrected(residualDataReal, correctedResiduals, file, extension="_model1")

            (epochs, yaw, errors) = estimation.solveLSEModel2(residualDataReal, orbitData, stationsData)
            correctedResiduals = utils.correctResData(orbitData, stationsData, residualDataReal, epochs, yaw)
            graphics.plotResidualsCorrected(residualDataReal, correctedResiduals, file, extension="_model2")

            (epochs, yaw, errors) = estimation.solveLSEModel3(residualDataReal, orbitData, stationsData)
            correctedResiduals = utils.correctResData(orbitData, stationsData, residualDataReal, epochs, yaw)
            graphics.plotResidualsCorrected(residualDataReal, correctedResiduals, file, extension="_model3")

            (epochs, yaw, errors) = estimation.solveLSEModel4(residualDataReal, orbitData, stationsData, Ne=Ne)
            correctedResiduals = utils.correctResData(orbitData, stationsData, residualDataReal, epochs, yaw)
            graphics.plotResidualsCorrected(residualDataReal, correctedResiduals, file, extension="_model4")
            '''
        except:
            print("Error on file {}".format(file))
