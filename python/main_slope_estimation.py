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
        graphics.plotSlopeOverBeta(slopeData, extension=file[0:3])
        graphics.plotInvertedOverBeta(slopeData, extension=file[0:3])
        graphics.plotEstimatedSlopeOverBeta(slopeData, extension=file[0:3], percentile=67)
        graphics.plotEstimatedSlopeOverBeta(slopeData, extension=file[0:3], percentile=95)
        graphics.plotEstimatedSlopeOverBeta(slopeData, extension=file[0:3], percentile=100)

