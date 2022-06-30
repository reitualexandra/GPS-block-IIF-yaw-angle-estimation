import os

# path to directory where res and orb files are stored
OUT = os.path.join("..", "residuals")

# path to generated figures directory
FIGS = os.path.join("..", "figs")

# path to generated residual directory
SIMS = os.path.join("..", "sim")

# path to GPSYAW.crd file - which contains list of stations used for residual generation
STATION_COORDINATE_LIST = os.path.join("..", "residuals", "GPSYAW.crd")

# elevation mask for measurements and other stuff
ELEVATION_MASK = 10
EARTH_RADIUS = 6371 * 1000

# GPS block IIF PCOs - in meters
PCO_x = 0.393
PCO_y = -0.017
PCO_z = 1.274