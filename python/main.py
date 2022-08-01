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
        print(file)
        try:
            file = str(file)
            year = file[3:5]
            doy = file[5:9]
            prn = file[11:14]
            man = file[9:11]
            man_simulated = "{}_sim".format(man)

            # Get stations data, orbit data, stations list and residuals for considered maneuver, doy, and prn
            orbfile = "YAW" + str(file[3:14]) + ".orb"  #YAW203580M1G25.orb
            orbitData = utils.getOrbData(orbfile)
            stationsList = utils.getStationsList(orbitData, stationsData)
            residualDataReal = utils.getResData(file)
            #residualDataReal = estimation.filterResiduals(residualDataReal)
            residualDataReal = utils.cleanResData(residualDataReal)
            Ne = len(orbitData['mjd']) - 1

            # Estimate yaw from real residuals
            graphics.plotResiduals(residualDataReal, stationsData, file)
            (epochs, yaw, errors) = estimation.solveLSEModel1(residualDataReal, orbitData, stationsData)
            #graphics.plotEstimatedYawInterpolated(epochs, yaw, orbitData, file, extension="_real_model1")
            graphics.plotEstimatedYawErrorbars(epochs, yaw, errors, orbitData, file,
                                               extension="_real_model1_bars")

            (epochs, yaw, errors) = estimation.solveLSEModel2(residualDataReal, orbitData, stationsData)
            #graphics.plotEstimatedYawInterpolated(epochs, yaw, orbitData, file, extension="_real_model2")
            graphics.plotEstimatedYawErrorbars(epochs, yaw, errors, orbitData, file,
                                               extension="_real_model2_bars")

            (epochs, yaw, errors) = estimation.solveLSEModel3(residualDataReal, orbitData, stationsData)
            #graphics.plotEstimatedYawInterpolated(epochs, yaw, orbitData, file, extension="_real_model3")
            graphics.plotEstimatedYawErrorbars(epochs, yaw, errors, orbitData, file,
                                               extension="_real_model3_bars")

            (epochs, yaw, errors) = estimation.solveLSEModel4(residualDataReal, orbitData, stationsData, Ne=Ne)
            #graphics.plotEstimatedYawInterpolated(epochs, yaw, orbitData, file, extension="_real_model4")
            graphics.plotEstimatedYawErrorbars(epochs, yaw, errors, orbitData, file,
                                               extension="_real_model4_bars")



        except:
            print("Error on file {}".format(file))
