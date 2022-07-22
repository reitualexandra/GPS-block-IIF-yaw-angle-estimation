import utils
import graphics
import estimation
import constants
import os

if __name__ == '__main__':

    files = os.listdir(constants.OUT)
    slope_files = [x for x in files if ("SlopeEst" in x)]

    for file in slope_files:
        slopeData = utils.getSlopeData(file)
        graphics.plotSlopeOverBeta(slopeData)
        graphics.plotInvertedOverBeta(slopeData)
        graphics.plotEstimatedSlopeOverBeta(slopeData)

