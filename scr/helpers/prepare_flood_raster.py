# Functions to help with preparing raster surface water maps, for example created by modified DeVries algorithm or downloaded from Copernicus GFM

import rasterio as rio
from rasterio import merge
from rasterio.mask import mask
from rasterio.merge import merge
import rioxarray
import xarray
import copy
import numpy as np
import pandas as pd
from datetime import datetime
from itertools import chain
import os

# Import module from prep_raster containing helpful functions to use
import importlib
import Helpers.visualize_modifiedDevries as visualize_DV
importlib.reload(visualize_DV)


# Function to mosaic raster tiles on a single date into one mosaic
def mosaic_rasters(tiles):
    
    # Open the raster and append to the list of tiles to mosaic
    raster_to_mosaic = []
    for i, tile in enumerate(tiles):
        raster = rio.open(tile)
        raster_to_mosaic.append(raster)
        if i!=(len(tiles)-1):
            raster = None
        
    print('Created list of tiles to mosaic')    

    # Use merge to create mosaic
    mosaic, output = merge(raster_to_mosaic)

    # Take a copy of and update the mosaic metadata
    output_meta = raster.meta.copy()
    output_meta.update(
        {"driver": "GTiff",
            "height": mosaic.shape[1],
            "width": mosaic.shape[2],
            "transform": output,
        }
    )

    return mosaic, output_meta
    
    
# Function for printing useful raster information
def print_raster(raster):
    print(
        f"shape: {raster.rio.shape}\n"
        f"resolution: {raster.rio.resolution()}\n"
        f"bounds: {raster.rio.bounds()}\n"
        f"sum: {raster.sum().item()}\n"
        f"CRS: {raster.rio.crs}\n"
        f"values: {np.unique(raster.values)}\n")

    
# Function to return geodataframe from bbox object
def bbox_to_gdf(lon0, lat0, lon1, lat1, CRS):
    poly = Polygon([[lon0, lat0],
                    [lon1,lat0],
                    [lon1,lat1],
                    [lon0, lat1]])
    gdf = gpd.GeoDataFrame(pd.DataFrame(['bbox'], columns = ['geom']), crs = CRS, geometry = [poly])
    return gdf


# Function to load raster using rioxarray into standard format with float datatype, 0 = background, 1 = water, no_data_value
def load_raster_for_comparison(raster_file, no_data_val):
    raster_orig = rioxarray.open_rasterio(raster_file) 
    raster = copy.copy(raster_orig.astype('float64'))
    raster.values[np.logical_and(raster.values!=0, raster.values!=1)]=no_data_val
    #print('Loaded raster unique values:', np.unique(raster.values))
    return raster

# Function to load raster using rioxarray into standard format with float datatype, 0 = background, 1 = water, no_data_value
def load_label_for_comparison(raster_file):
    raster_orig = rioxarray.open_rasterio(raster_file) 
    raster = copy.copy(raster_orig.astype('float64'))
    raster.values[np.logical_or(raster.values==1, raster.values==2)]=1
    raster.values[raster.values!=1]=0
    #print('Loaded raster unique values:', np.unique(raster.values))
    return raster

# Function to reproject match one raster to the extent and resolution of another raster, including assigning the coordinates
def reproj_match_raster(raster_to_reproj, raster_to_match):
    raster_repr_match = raster_to_reproj.rio.reproject_match(raster_to_match)
    raster_repr_match = raster_repr_match.assign_coords({'x': raster_to_match.x, 'y': raster_to_match.y,})
    return raster_repr_match

# Function to get the Devries map to use as bbox for plotting other images with larger extent
def get_Devries_rasterBbox(Devries_mosaic, targetCRS):
    raster_orig = rioxarray.open_rasterio(Devries_mosaic) 
    raster_orig = raster_orig.rio.reproject(targetCRS)
    return raster_orig

# Function to get the bbox in the rasterio.bbox format for a raster loaded in using rioxarray
def get_bbox_rioxarray(raster_rioxarray):
    boundingbox = raster_rioxarray.rio.bounds()
    bbox = rio.coords.BoundingBox(boundingbox[0], boundingbox[1], boundingbox[2], boundingbox[3])
    return bbox


# Function to get the matching DeVries and GFM raster for a given date, with the GFM raster being reprojected and matched to the DeVries map
def get_Devries_and_GFM_rasters(date, Devries_mosaics, GFM_mosaics):

    # Get the image dates
    GFM_dates = [str(datetime.strptime(mosaic.stem[-10:].replace('_',''), '%Y%m%d').date()) for mosaic in GFM_mosaics]
    Devries_dates = [str(datetime.strptime(mosaic.stem[-14:-6], '%Y%m%d').date()) for mosaic in Devries_mosaics]
    
    # Check if date exists for both images
    if date in GFM_dates and date in Devries_dates:
        print('Images exist for target date')
    elif date in GFM_dates and date not in Devries_dates:
        print('Modified DeVries image does not exist for target date. STOPPING')
        return 0, 0
    elif date not in GFM_dates and date in Devries_dates:
        print('GFM image does not exist for target date. STOPPING')
        return 0, 0
    else:
        print('Images do not exist for target date. STOPPING')
        return 0, 0
        
    # Load the rasters in float datatype, with 0 (no water), 1 (water) and -1 (no values)
    raster_GFM = load_raster_for_comparison(GFM_mosaics[GFM_dates.index(date)], -10)
    raster_Devries = load_raster_for_comparison(Devries_mosaics[Devries_dates.index(date)], -20)

    # Show information about the rasters
    # print('Raster Devries')
    # print_raster(raster_Devries)
    # print('Original Raster GFM')
    # print_raster(raster_GFM)

    # Reproject match the GFM raster to the extent and resolution of the Devries raster
    raster_GFM_repr_match = reproj_match_raster(raster_GFM, raster_Devries)

    # Show the information about the reprojected raster (GFM) and check it matches the target raster (DeVries)
    # print('Reprojected and matched Raster GFM')
    # print_raster(raster_GFM_repr_match)
    
    return raster_Devries, raster_GFM_repr_match

# Function to get the UNOSAT raster, reprojected and matched to the DeVries map
def get_UNOSAT_raster(date, UNOSAT_mosaics, Devries_mosaic):

    # Get the image dates
    UNOSAT_dates = ['20220619','20220525']
    
    # Check if date exists for both images
    if date in UNOSAT_dates:
        print('UNOSAT image exists for target date')
    else:
        print('UNOSAT image does not exist for target date. STOPPING')
        return 0, 0
        
    # Load the rasters in float datatype, with 0 (no water), 1 (water) and -1 (no values)
    raster_UNOSAT = load_raster_for_comparison(UNOSAT_mosaics[UNOSAT_dates.index(date)], -10)
    raster_Devries = load_raster_for_comparison(Devries_mosaic, -20)

    # Show information about the rasters
    # print('Raster Devries')
    # print_raster(raster_Devries)
    # print('Original Raster GFM')
    # print_raster(raster_GFM)

    # Reproject match the GFM raster to the extent and resolution of the Devries raster
    raster_UNOSAT_repr_match = reproj_match_raster(raster_UNOSAT, raster_Devries)

    # Show the information about the reprojected raster (GFM) and check it matches the target raster (DeVries)
    # print('Reprojected and matched Raster GFM')
    # print_raster(raster_GFM_repr_match)
    
    return raster_UNOSAT_repr_match


# Function to take rioxarray and convert to numpy array for computing confusion matrix, sets all non water values to 0, makes bindary integer
def raster_to_binary_list(raster):
    array = raster.astype('int').to_numpy()
    # array[array!=1]=0
    array[np.logical_and(array!=1,array!=0)]=255
    array_list = list(chain.from_iterable(list(chain.from_iterable(array.tolist()))))
    return array_list

# Function to convert raster from original to target CRS and write as a tif in the same folder
# raster should be the  raster file including complete file path
def write_reproj_raster(raster, target_crs, suffix_to_add):

    with rio.open(raster) as src:
        transform, width, height = calculate_default_transform(
            src.crs, target_crs, src.width, src.height, *src.bounds)
        kwargs = src.meta.copy()
        kwargs.update({
            'crs': target_crs,
            'transform': transform,
            'width': width,
            'height': height
        })

        with rio.open(str(raster).replace('.tif', suffix_to_add+'.tif'), 'w', **kwargs) as dst:
            for i in range(1, src.count + 1):
                reproject(
                    source=rio.band(src, i),
                    destination=rio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=dst_crs,
                    resampling=Resampling.nearest)
                
# Function to get the matching DeVries and GFM raster for a given date, with the GFM raster being reprojected and matched to the DeVries map
def get_Devries_GFM_IMPACT_rasters(date, Devries_mosaics, GFM_mosaics, IMPACT_mosaics):

    # Get the image dates
    GFM_dates = [str(datetime.strptime(mosaic.stem[-10:].replace('_',''), '%Y%m%d').date()) for mosaic in GFM_mosaics]
    Devries_dates = [str(datetime.strptime(mosaic.stem[-14:-6], '%Y%m%d').date()) for mosaic in Devries_mosaics]
    IMPACT_dates = [str(datetime.strptime(file.stem, '%Y%m%d').date()) for file in IMPACT_mosaics]
    
    # Check if date exists for both images
    if date in GFM_dates and date in Devries_dates:
        print('Images exist for target date')
    elif date in GFM_dates and date not in Devries_dates:
        print('Modified DeVries image does not exist for target date. STOPPING')
        return 0, 0
    elif date not in GFM_dates and date in Devries_dates:
        print('GFM image does not exist for target date. STOPPING')
        return 0, 0
    else:
        print('Images do not exist for target date. STOPPING')
        return 0, 0
        
    # Load the rasters in float datatype, with 0 (no water), 1 (water) and -1 (no values)
    raster_GFM = load_raster_for_comparison(GFM_mosaics[GFM_dates.index(date)], -10)
    raster_Devries = load_raster_for_comparison(Devries_mosaics[Devries_dates.index(date)], -20)
    raster_IMPACT = load_raster_for_comparison(IMPACT_mosaics[IMPACT_dates.index(date)], -30)

    # Show information about the rasters
    # print('Raster Devries')
    # print_raster(raster_Devries)
    # print('Original Raster GFM')
    # print_raster(raster_GFM)

    # Reproject match the GFM raster to the extent and resolution of the Devries raster
    raster_GFM_repr_match = reproj_match_raster(raster_GFM, raster_Devries)
    raster_IMPACT_repr_match = reproj_match_raster(raster_IMPACT, raster_Devries)

    # Show the information about the reprojected raster (GFM) and check it matches the target raster (DeVries)
    # print('Reprojected and matched Raster GFM')
    # print_raster(raster_GFM_repr_match)
    
    return raster_Devries, raster_GFM_repr_match, raster_IMPACT_repr_match

# Function to resample and match (align) one raster to another target raster, writes raster out as tif
# target=the one we want to match to
# source=the one to modify to match to the source raster
def matchWriteRaster(sourceRaster, targetRaster, outputFile):
    with rioxarray.open_rasterio(sourceRaster) as sourceRaster_r:
        with rioxarray.open_rasterio(targetRaster) as targetRaster_r:

            sourceRaster_crs = str(sourceRaster_r.rio.crs)
            targetRaster_gt = targetRaster_r.rio.transform()
            targetRaster_res = (targetRaster_gt[0], -targetRaster_gt[4])
            targetRaster_crs = str(targetRaster_r.rio.crs)
            
            raster_extent = targetRaster_r.rio.bounds()
            command = f"gdalwarp \
                                -s_srs {sourceRaster_crs} \
                                -t_srs {targetRaster_crs} \
                                -dstnodata {str(targetRaster_r.rio.nodata)} \
                                -tr {str(targetRaster_res[0])} {str(targetRaster_res[1])} \
                                -r average \
                                -te {str(raster_extent[0])} {str(raster_extent[1])} {str(raster_extent[2])} {str(raster_extent[3])} \
                                -te_srs {targetRaster_crs} \
                                -ot Float64 \
                                -of GTiff {str(sourceRaster)} \
                                {str(outputFile)}"
            os.system(command)
                                # -overwrite \

# Function to load raster using rioxarray into standard format with float datatype, FIA values in range 0 to 1
def load_MatchedRaster_for_comparison(raster_file, no_data_val):
    raster_orig = rioxarray.open_rasterio(raster_file) 
    raster = copy.copy(raster_orig.astype('float64'))
    raster.values[np.logical_or(raster.values<0, raster.values>1)]=no_data_val
    #print('Loaded raster unique values:', np.unique(raster.values))
    return raster


# Function to get the matching DeVries and GFM raster for a given date, with the GFM raster being reprojected and matched to the DeVries map
def get_Devries_GFM_IMPACT_MatchedRasters(date, Devries_mosaics, GFM_mosaics, IMPACT_mosaics, targetRaster):

    # Get the image dates
    GFM_dates = [str(datetime.strptime(mosaic.stem, '%Y%m%d').date()) for mosaic in GFM_mosaics]
    Devries_dates = [str(datetime.strptime(mosaic.stem, '%Y%m%d').date()) for mosaic in Devries_mosaics]
    IMPACT_dates = [str(datetime.strptime(file.stem, '%Y%m%d').date()) for file in IMPACT_mosaics]
    
    # Check if date exists for both images
    if date in GFM_dates and date in Devries_dates:
        print('Images exist for target date')
    elif date in GFM_dates and date not in Devries_dates:
        print('Modified DeVries image does not exist for target date. STOPPING')
        return 0, 0
    elif date not in GFM_dates and date in Devries_dates:
        print('GFM image does not exist for target date. STOPPING')
        return 0, 0
    else:
        print('Images do not exist for target date. STOPPING')
        return 0, 0
        
    # Load the rasters in float datatype, with 0 (no water), 1 (water) and -1 (no values)
    raster_GFM = load_MatchedRaster_for_comparison(GFM_mosaics[GFM_dates.index(date)], 255)
    raster_Devries = load_MatchedRaster_for_comparison(Devries_mosaics[Devries_dates.index(date)], 255)
    raster_IMPACT = load_MatchedRaster_for_comparison(IMPACT_mosaics[IMPACT_dates.index(date)], 255)
    
    raster = load_MatchedRaster_for_comparison(targetRaster, 255)

    # Reproject match the GFM raster to the extent and resolution of the Devries raster
    raster_GFM_repr_match = reproj_match_raster(raster_GFM, raster)
    raster_Devries_repr_match = reproj_match_raster(raster_Devries, raster)
    raster_IMPACT_repr_match = reproj_match_raster(raster_IMPACT, raster)

    return raster_Devries_repr_match, raster_GFM_repr_match, raster_IMPACT_repr_match


# Function to get the matching DeVries and GFM raster for a given date, with the GFM raster being reprojected and matched to the DeVries map
def get_Fusion_GFD_MatchedRasters(date, Fusion_mosaics, GFD_mosaics, targetRaster):

    # Get the image dates
    Fusion_dates = [str(datetime.strptime(mosaic.stem, '%Y%m%d').date()) for mosaic in Fusion_mosaics]
    GFD_dates = [str(datetime.strptime(mosaic.stem, '%Y%m%d').date()) for mosaic in GFD_mosaics]
    
    # Check if date exists for both images
    if date in Fusion_dates and date in GFD_dates:
        print('Images exist for target date')
    elif date in Fusion_dates and date not in GFD_dates:
        print('GFD image does not exist for target date. STOPPING')
        return 0, 0
    elif date not in Fusion_dates and date in GFD_dates:
        print('Fusion image does not exist for target date. STOPPING')
        return 0, 0
    else:
        print('Images do not exist for target date. STOPPING')
        return 0, 0
        
    # Load the rasters in float datatype, with 0 (no water), 1 (water) and -1 (no values)
    raster_Fusion = load_MatchedRaster_for_comparison(Fusion_mosaics[Fusion_dates.index(date)], 255)
    raster_GFD = load_MatchedRaster_for_comparison(GFD_mosaics[GFD_dates.index(date)], 255)

    raster = load_MatchedRaster_for_comparison(targetRaster, 255)
    
    # Reproject match (probably not required)
    raster_Fusion_repr_match = reproj_match_raster(raster_Fusion, raster)
    raster_GFD_repr_match = reproj_match_raster(raster_GFD, raster)

    return raster_Fusion_repr_match, raster_GFD_repr_match