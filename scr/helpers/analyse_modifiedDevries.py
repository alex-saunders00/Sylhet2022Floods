# Functions to help with analysing the outputs from the modified Devries algorithm surface water maps

from pathlib import Path
import os
import sys
import pandas as pd
import numpy as np
import rasterio
from rasterio.plot import show
from rasterio.mask import mask
import pycrs
from shapely.geometry import mapping
from itertools import chain
from datetime import datetime
import copy
import collections
from functools import partial
import pyproj
from shapely.ops import transform
from shapely.geometry import Point
from sklearn import metrics
import rioxarray


# Import module from prep_raster containing helpful functions to use
import importlib
import Helpers.prepare_flood_raster as prep_raster
importlib.reload(prep_raster)

# Create a function to create a raster of values 10, the same dimensions as the flood raster, to be used for clipping to the ROI later - ONLY NEEDS TO BE RUN ONCE
def create_temp_raster_for_mask(mosaic, rootPath, inputPath):
    
    # Load the flood map raster - simply takes the first one of all the dates, assuming they have same extent at each date
    raster_tmp = rasterio.open(mosaic)
    tmp_data = raster_tmp.read(1)

    # Set all values to 10 and reshape from 2D to 3D
    tmp_data[:]=10
    np.reshape(tmp_data, tmp_data.shape + (1,))

    # Write out converting from numpy array back to geotiff
    # convert the numpy array back to raster file
    with rasterio.open(rootPath/inputPath/'temp.tif', 'w',
                       driver='GTiff',
                       height=tmp_data.shape[0],
                       width=tmp_data.shape[1],
                       count=1,
                       dtype=tmp_data.dtype,
                       crs=raster_tmp.crs,
                       transform=raster_tmp.transform) as dst:
        dst.write(tmp_data, indexes=1)
        
# Create a function to create a raster of values 10, the same dimensions as the flood raster, to be used for clipping to the ROI later - ONLY NEEDS TO BE RUN ONCE
def create_temp_raster_for_mask_wCRS(mosaic, rootPath, inputPath, targetCRS):
    
    # Load the flood map raster - simply takes the first one of all the dates, assuming they have same extent at each date
    raster_tmp = rasterio.open(mosaic)
    tmp_data = raster_tmp.read(1)

    # Set all values to 10 and reshape from 2D to 3D
    tmp_data[:]=10
    np.reshape(tmp_data, tmp_data.shape + (1,))

    # Write out converting from numpy array back to geotiff
    # convert the numpy array back to raster file
    with rasterio.open(rootPath/inputPath/'temp.tif', 'w',
                       driver='GTiff',
                       height=tmp_data.shape[0],
                       width=tmp_data.shape[1],
                       count=1,
                       dtype=tmp_data.dtype,
                       crs=targetCRS,
                       transform=raster_tmp.transform) as dst:
        dst.write(tmp_data, indexes=1)

        
# Create function to create the raster mask for the ROI, using the temp raster containing all values 10
def create_roi_mask(raster_tmp, roi): 
    #roi_newCRS = roi.to_crs(crs=raster_tmp.crs.data)
    roi_newCRS = roi.to_crs(crs=int(raster_tmp.crs.data['init'][5:]))
    roi_mask, roi_mask_transform = mask(raster_tmp, roi_newCRS.geometry.apply(mapping), crop=False)
    out_meta = raster_tmp.meta.copy()
    epsg_code = int(raster_tmp.crs.data['init'][5:])
    out_meta.update({"driver": "GTiff",
                    "height": roi_mask.shape[1],
                    "width": roi_mask.shape[2],
                    "transform": roi_mask_transform,
                    "crs": pycrs.parse.from_epsg_code(epsg_code).to_proj4()}) # replace +init=epsg:4326 with EPSG:4326 
    roi_mask = roi_mask[0] # reduce to 2D
    return roi_mask


# Create a function to compute the fractional area inundated for a given date
def get_fracArea(mosaic, roi_mask):
    
    # Get the dates as date format
    mosaic_date = str(datetime.strptime(mosaic.stem[-14:-6], '%Y%m%d').date())

    # Load the raster
    with rasterio.open(mosaic) as file:
        raster = file.read(1)

    # Set the water values to 1
    water = copy.copy(raster)
    water[water==1]=1

    # # Get the number of pixels of water for the entire image
    # water_vals = list(chain.from_iterable(water.tolist()))
    # counter = collections.Counter(water_vals) # 0 = none, 1 = water, other = invalid
    # pctWater = counter[1] / (counter[0] + counter[1]) #sum(counter.values()
    # #print(counter)
    # print('Image % water at {0}: {1:2.1f}%'.format(mosaic_date, pctWater*100))
    
    # Apply the mask to the water_sum raster by performing an array sum operation
    water_roi = np.add(water, roi_mask)

    # Get the number of pixels of water for each date for the ROI only (pixel values in ROI have been shifted by +10 from original values)
    water_roi_vals = list(chain.from_iterable(water_roi.tolist()))
    counter_roi = collections.Counter(water_roi_vals) # 10 = none, 11 = water, other = invalid or outside ROI
    pctWater_roi = counter_roi[11] / (counter_roi[10] + counter_roi[11])
    #print(counter_roi)
    print('ROI % water at {0}: {1:2.1f}%'.format(mosaic_date, pctWater_roi*100))
    
    return mosaic_date, pctWater_roi

# Create a function to get the image and ROI % water for all image dates and a given ROI
def get_fracArea_timeSeries(mosaics, raster_tmp, roi):

    # Run the function to create the raster mask for the ROI, using the temp raster containing all values 10
    roi_mask = create_roi_mask(raster_tmp, roi)

    # Loop over the image dates and run the function to compute the fractional area inundated for a provided raster on a given date
    mosaic_dates = []
    pctWaters = []

    for n, mosaic in enumerate(mosaics):

        # Get the fractional area inundated
        mosaic_date, pctWater_roi = get_fracArea(mosaic, roi_mask)

        # Record the outputs for the given date
        mosaic_dates.append(mosaic_date)
        pctWaters.append(pctWater_roi)
        
    # Create dataframe of the results - fractional flooded area
    fracFloodArea = pd.DataFrame(data = [mosaic_dates, pctWaters], index = ['Date','FFA']).T
    fracFloodArea.set_index('Date', inplace=True)
    fracFloodArea.index = pd.to_datetime(fracFloodArea.index)

    return fracFloodArea


# Create a function to compute the fractional area inundated for a given date, for the GFM maps
def get_fracArea_GFM(mosaic, roi_mask, raster_to_match):
    
    # Get the dates as date format
    mosaic_date = str(datetime.strptime(mosaic.stem[-10:].replace('_',''), '%Y%m%d').date())

    # Load the raster using rioxarray
    raster = prep_raster.load_raster_for_comparison(mosaic, -10)
    
    # Reproject and match the GFM raster to the modified Devries raster to make the processing the same
    raster_reproj = prep_raster.reproj_match_raster(raster, raster_to_match)

    # Set the water values to 1
    water = copy.copy(raster_reproj)
    water.values[water.values==1]=1
    water.values[water.values!=1]=0
  
    # Apply the mask to the water_sum raster by performing an array sum operation
    water_roi = np.add(water, roi_mask)

    # Get the number of pixels of water for each date for the ROI only (pixel values in ROI have been shifted by +10 from original values)
    water_roi_vals = list(chain.from_iterable(list(chain.from_iterable(water_roi.values.tolist()))))
    counter_roi = collections.Counter(water_roi_vals) # 10 = none, 11 = water, other = invalid or outside ROI
    pctWater_roi = counter_roi[11] / (counter_roi[10] + counter_roi[11])
    #print(counter_roi)
    print('ROI % water at {0}: {1:2.1f}%'.format(mosaic_date, pctWater_roi*100))
    
    return mosaic_date, pctWater_roi

# Create a function to get the image and ROI % water for all image dates and a given ROI, for the Copernicus GFM data
def get_fracArea_timeSeries_GFM(mosaics, raster_tmp, raster_to_match, roi):

    # Run the function to create the raster mask for the ROI, using the temp raster containing all values 10
    roi_mask = create_roi_mask(raster_tmp, roi)

    # Loop over the image dates and run the function to compute the fractional area inundated for a provided raster on a given date
    mosaic_dates = []
    pctWaters = []

    for n, mosaic in enumerate(mosaics):

        # Get the fractional area inundated
        mosaic_date, pctWater_roi = get_fracArea_GFM(mosaic, roi_mask, raster_to_match)

        # Record the outputs for the given date
        mosaic_dates.append(mosaic_date)
        pctWaters.append(pctWater_roi)
        
    # Create dataframe of the results - fractional flooded area
    fracFloodArea = pd.DataFrame(data = [mosaic_dates, pctWaters], index = ['Date','FFA']).T
    fracFloodArea.set_index('Date', inplace=True)
    fracFloodArea.index = pd.to_datetime(fracFloodArea.index)

    return fracFloodArea


# Create function to perform geodesic point buffer, using lat/lon input and a transformation to an azimuthal equidistant projection, which has all points at proportionally correct distances from the center point
def geodesic_point_buffer(lat, lon, km):
    # WGS84 and azimuthal equidistant projection
    proj_wgs84 = pyproj.Proj('+proj=longlat +datum=WGS84')
    aeqd_proj = '+proj=aeqd +lat_0={lat} +lon_0={lon} +x_0=0 +y_0=0'
    project = partial(
        pyproj.transform,
        pyproj.Proj(aeqd_proj.format(lat=lat, lon=lon)),
        proj_wgs84)
    buf = Point(0, 0).buffer(km * 1000)  # distance in metres
    return transform(project, buf).exterior.coords[:]


# Function to get confusion matrix and other metrics for comparison of two water rasters
def get_confusion_matrix_metrics(raster_list_actual, raster_list_prediction):
    
    # Compute the confusion matrix
    # TP: DeVries and GFM = water, TN: DeVries and GFM = no water, FP: Only GFM = water, FN: Only DeVries = water
    tn, fp, fn, tp = metrics.confusion_matrix(raster_list_actual, raster_list_prediction).ravel()

    # Compute precision, recall, accuracy, intersection over union
    pre = tp / (tp+fp)
    rec = tp / (tp+fn)
    acc = (tp+tn) / (tp+fp+fn+tn)
    iou = tp / (tp+fp+fn) 
    
    metrics_list = [tn, fp, fn, tp, pre, rec, acc, iou]
    
    return metrics_list


# Function to get list of FIA values from the reprojected and matched raster at 500 m resolution, clipped within the ROI
# Raster should be already opened with rioxarray
def getFIAValsROI(rasterFile, roi, rasterBand=0):
    
    # Open the ratser with rioxarray and get first band
    raster = rioxarray.open_rasterio(rasterFile)[rasterBand]
    
    # Clip raster to ROI geometry
    clipRaster = raster.rio.clip([roi]) # roi is a shapely Polygon
    
    # Get a list of all valid pixel values, ignoring / remove nodata values
    FIAVals = clipRaster.values[np.logical_and(clipRaster.values>=0, clipRaster.values<=1)]
    
    return FIAVals

# Function to get stats (mean, median, stdev, 5 and 95%iles) from FIA values list
def getFIAStats(FIAVals):
    FIAStats =[np.mean(FIAVals),
    np.median(FIAVals),
    np.std(FIAVals),
    np.quantile(FIAVals, 0.05),
    np.quantile(FIAVals, 0.95)]
    return FIAStats