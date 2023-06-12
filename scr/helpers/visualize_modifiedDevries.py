# Functions to help with visualizing the outputs from the modified Devries algorithm surface water maps

from pathlib import Path
import os
import sys
import pandas as pd
import numpy as np
import rasterio
from itertools import chain
from datetime import datetime
import math
from pyproj import Transformer
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import copy
import collections
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patheffects as pe
import matplotlib.font_manager as fm
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import scipy

# Import module from prep_raster containing helpful functions to use
import importlib
import Helpers.prepare_flood_raster as prep_raster
importlib.reload(prep_raster)


# Functions for rounding to get neat coordinate labels
def round_nearest(num: float, to: float) -> float:
    return round(num / to) * to  # Credited to Paul H.

def round_down(num: float, to: float) -> float:
    nearest = round_nearest(num, to)
    if math.isclose(num, nearest): return num
    return nearest if nearest < num else nearest - to

def round_up(num: float, to: float) -> float:
    nearest = round_nearest(num, to)
    if math.isclose(num, nearest): return num
    return nearest if nearest > num else nearest + to

        
# Create function to get ticks and labels for coordinate for raster maps
def create_coord_ticks(bbox, orig_crs, labels_crs, rounding_x, rounding_y): 
# bbox is the bounds of the raster, orig_crs is the original CRS, labels_crs is a different target CRS for the coord labels if desired
# rounding_x / rounding_y is the rounding multiple in the labels_crs for x / y
# Specify CRS in this format 'epsg:32646', 'epsg:4326'

    # Get the min and max extent and transform to latlon EPSG4326
    transformer_to = Transformer.from_crs(orig_crs, labels_crs) # from, to
    x_latlon, y_latlon = transformer_to.transform([bbox.left, bbox.right], [bbox.bottom, bbox.top])

    # Create the tick labels as coordinates within the extent, rounded to the nearest 0.1 degrees
    x_latlon1 = round_up(min(x_latlon), rounding_x)
    x_latlonn = round_up(max(x_latlon), rounding_x)
    xticks_latlon = np.arange(x_latlon1, x_latlonn + rounding_x, rounding_x).tolist()
    y_latlon1 = round_up(min(y_latlon), rounding_y)
    y_latlonn = round_up(max(y_latlon), rounding_y)
    yticks_latlon = np.arange(y_latlon1, y_latlonn + rounding_y, rounding_y).tolist()

    # x and y must be of same length for converting back to original CRS, so pad out with a repeat value for whichever is shorter of x and y
    x_len_orig = len(xticks_latlon)
    y_len_orig = len(yticks_latlon)
    
    if len(xticks_latlon) < len(yticks_latlon):
        for i in range(0, (len(yticks_latlon) - len(xticks_latlon))):
            xticks_latlon.append(xticks_latlon[len(xticks_latlon)-1])
    elif len(xticks_latlon) > len(yticks_latlon):
        for i in range(0, (len(xticks_latlon) - len(yticks_latlon))):
            yticks_latlon.append(yticks_latlon[len(yticks_latlon)-1])

    # Convert the latlon tick labels back to ticks in the original CRS
    transformer_back = Transformer.from_crs(labels_crs, orig_crs) # from, to
    xticks_origCRS, yticks_origCRS = transformer_back.transform(xticks_latlon, yticks_latlon)
    
    # Remove duplicates
    # xticks_origCRS = [*set(xticks_origCRS)]
    # xticks_latlon = [*set(xticks_latlon)]
    # yticks_origCRS = [*set(yticks_origCRS)]
    # yticks_latlon = [*set(yticks_latlon)]
    
    # # Return to the original length, to avoid duplicate labels
    # xticks_origCRS = xticks_origCRS[:x_len_orig-1]
    # xticks_latlon = xticks_latlon[:y_len_orig-1]
    # yticks_origCRS = yticks_origCRS[:y_len_orig-1]
    # yticks_latlon = yticks_latlon[:x_len_orig-1]
        
    # Return the needed outputs
    return xticks_origCRS, yticks_origCRS, xticks_latlon, yticks_latlon


# Create function to get ticks and labels for coordinate for raster maps
def create_coord_ticks_region(bbox, orig_crs, labels_crs, rounding_x, rounding_y): 
# bbox is the bounds of the raster, orig_crs is the original CRS, labels_crs is a different target CRS for the coord labels if desired
# rounding_x / rounding_y is the rounding multiple in the labels_crs for x / y
# Specify CRS in this format 'epsg:32646', 'epsg:4326'

    # Get the min and max extent and transform to latlon EPSG4326
    transformer_to = Transformer.from_crs(orig_crs, labels_crs) # from, to
    x_latlon, y_latlon = transformer_to.transform([bbox.left[0], bbox.right[0]], [bbox.bottom[0], bbox.top[0]])

    # Create the tick labels as coordinates within the extent, rounded to the nearest 0.1 degrees
    x_latlon1 = round_up(min(x_latlon), rounding_x)
    x_latlonn = round_up(max(x_latlon), rounding_x)
    xticks_latlon = np.arange(x_latlon1, x_latlonn + rounding_x, rounding_x).tolist()
    y_latlon1 = round_up(min(y_latlon), rounding_y)
    y_latlonn = round_up(max(y_latlon), rounding_y)
    yticks_latlon = np.arange(y_latlon1, y_latlonn + rounding_y, rounding_y).tolist()

    # x and y must be of same length for converting back to original CRS, so pad out with a repeat value for whichever is shorter of x and y
    x_len_orig = len(xticks_latlon)
    y_len_orig = len(yticks_latlon)
    
    if len(xticks_latlon) < len(yticks_latlon):
        for i in range(0, (len(yticks_latlon) - len(xticks_latlon))):
            xticks_latlon.append(xticks_latlon[len(xticks_latlon)-1])
    elif len(xticks_latlon) > len(yticks_latlon):
        for i in range(0, (len(xticks_latlon) - len(yticks_latlon))):
            yticks_latlon.append(yticks_latlon[len(yticks_latlon)-1])

    # Convert the latlon tick labels back to ticks in the original CRS
    transformer_back = Transformer.from_crs(labels_crs, orig_crs) # from, to
    xticks_origCRS, yticks_origCRS = transformer_back.transform(xticks_latlon, yticks_latlon)
    
    # Remove duplicates
    # xticks_origCRS = [*set(xticks_origCRS)]
    # xticks_latlon = [*set(xticks_latlon)]
    # yticks_origCRS = [*set(yticks_origCRS)]
    # yticks_latlon = [*set(yticks_latlon)]
    
    # # Return to the original length, to avoid duplicate labels
    # xticks_origCRS = xticks_origCRS[:x_len_orig-1]
    # xticks_latlon = xticks_latlon[:y_len_orig-1]
    # yticks_origCRS = yticks_origCRS[:y_len_orig-1]
    # yticks_latlon = yticks_latlon[:x_len_orig-1]
        
    # Return the needed outputs
    return xticks_origCRS, yticks_origCRS, xticks_latlon, yticks_latlon


# Function to plot simple Devries (or GFM) map, for doing comparison plots of Devries vs GFM
def plot_simple_raster_for_comparison(raster, ax, cmap, norm, title, orig_crs, labels_crs, rounding_x, rounding_y):
    
    bbox = prep_raster.get_bbox_rioxarray(raster)
    raster_plot = ax.imshow(np.squeeze(raster), cmap=cmap, norm=norm, 
                            extent = [bbox.left, bbox.right, bbox.bottom, bbox.top], 
                            interpolation='none')
    ax.set_title(title, fontsize=10)
    xticks_origCRS, yticks_origCRS, xticks_latlon, yticks_latlon = create_coord_ticks(bbox, orig_crs, labels_crs, rounding_x, rounding_y)
    # ax.xticks(xticks_origCRS, ['{:.1f}'.format(x) +'\N{DEGREE SIGN}E' for x in yticks_latlon], fontsize=6)
    # ax.yticks(yticks_origCRS, ['{:.1f}'.format(x) +'\N{DEGREE SIGN}N' for x in xticks_latlon], fontsize=6)    
    ax.set_xticks(xticks_origCRS)
    ax.set_xticklabels(['{:.1f}'.format(x) +'\N{DEGREE SIGN}E' for x in yticks_latlon], fontsize=6)
    ax.set_yticks(yticks_origCRS)
    ax.set_yticklabels(['{:.1f}'.format(x) +'\N{DEGREE SIGN}N' for x in xticks_latlon], fontsize=6)   
    ax.set_xlim(bbox.left, bbox.right)
    ax.set_ylim(bbox.bottom, bbox.top)

    # Add legend
    legend_labels = {'xkcd:azure': 'Water',
                     'lightgray': 'Background', 
                     'black': 'No values'}
    patches = [Patch(color=color, label=label) for color, label in legend_labels.items()]
    ax.legend(handles=patches, #bbox_to_anchor=(1.35, 1),
              facecolor='white', fontsize=8)

    
    
# Function to plot the water across two dates, performs a raster sum to identify agreement / disagreement in water extent
def map_water_diff_dates_subplot(source, raster_1, raster_2, date1, date2, orig_crs, labels_crs, rounding_x, rounding_y, ax, cmap, cols_to_use, norm, raster_mask, districts, upazillas):
    
    # Set the water values in date1 to 1, and in date2 to 2
    water_1 = copy.deepcopy(raster_1)
    water_1.values[water_1.values<0]=0
    water_1.values[water_1.values==1]=1

    water_2 = copy.deepcopy(raster_2)
    water_2.values[water_2.values<0]=0
    water_2.values[water_2.values==1]=2
    
    # Sum the water values from the two rasters
    water_sum = np.add(water_1, water_2)

#     # Apply the mask to the water_sum raster by performing an array sum operation
#     water_sum_roi = np.add(water_sum, raster_mask)
#     np.unique(water_sum_roi)

#     # Get the number of pixels of water for each date for the ROI only (pixel values in ROI have been shifted by +10 from original values)
#     water_sum_roi_vals = list(chain.from_iterable(list(chain.from_iterable(water_sum_roi.values.tolist()))))
#     counter_roi = collections.Counter(water_sum_roi_vals) # 10 = none, 11 = date1 only, 12 = date2 only, 13 = date1 and date2
#     pctWater_roi_1 = (counter_roi[11] + counter_roi[13]) / (counter_roi[10] + counter_roi[11] + counter_roi[12] + counter_roi[13])
#     pctWater_roi_2 = (counter_roi[12] + counter_roi[13]) / (counter_roi[10] + counter_roi[11] + counter_roi[12] + counter_roi[13])
#     # print(counter_roi)
#     # print('ROI % water at {0}: {1:2.1f}%'.format(date1, pctWater_roi_1*100))
#     # print('ROI % water at {0}: {1:2.1f}%'.format(date2, pctWater_roi_2*100))

    # Create plot
    bbox = prep_raster.get_bbox_rioxarray(water_sum)
    raster_plot = ax.imshow(np.squeeze(water_sum), cmap=cmap, norm=norm, 
                            extent = [bbox.left, bbox.right, bbox.bottom, bbox.top], interpolation='none', zorder=1)

    # Add the Sylhet District outline
    #sylhet_dist.to_crs(epsg=32646).boundary.plot(ax=ax, color = 'yellow', linewidth=2, label='Sylhet District', zorder=3)
    districts.to_crs(epsg=32646).boundary.plot(ax=ax, color = 'yellow', linewidth=1, zorder=2)
    #sylhet_upa.to_crs(epsg=32646).plot(ax=ax, color = 'yellow', alpha=0.3, linewidth=1, label='Sylhet Upazilla', zorder=3)
    upazillas.to_crs(epsg=32646).boundary.plot(ax=ax, color = 'yellow', linewidth=0.5, zorder=2)

    # Set title and legend
    #ax.set_title(str(source)+': '+str(date1)+' vs '+str(date2))
    ax.set_title(str(source), fontsize=20)
    # legend_labels = {cols_to_use[1] : 'Background',
    #                  cols_to_use[2]: str(date1)+' only',
    #                  cols_to_use[3]: str(date2)+' only',
    #                  cols_to_use[4]: 'Both dates'}
    # patches = [Patch(color=color, label=label) for color, label in legend_labels.items()]
    # ax.legend(handles=patches, #bbox_to_anchor=(1.35, 1),
    #           facecolor='white')

    # Set the ticks and labels
    # Run the function to get the ticks and ticklabels
    xticks_origCRS, yticks_origCRS, xticks_latlon, yticks_latlon = create_coord_ticks(bbox, orig_crs, labels_crs, rounding_x, rounding_y)
    ax.set_xticks(xticks_origCRS, ['{:.1f}'.format(x) +'\N{DEGREE SIGN}E' for x in yticks_latlon], fontsize=12)
    ax.set_yticks(yticks_origCRS, ['{:.1f}'.format(x) +'\N{DEGREE SIGN}N' for x in xticks_latlon], fontsize=12)
    ax.set_xlim(bbox.left, bbox.right)
    ax.set_ylim(bbox.bottom, bbox.top)
    # ax.tick_params(top=False, bottom=False, left=False, right=False,
    #             labelleft=False, labelbottom=False)
    
#     # Add text box with the fractional area inundated
#     textstr = '\n'.join(('ROI % water',
#         ('{0}: {1:2.1f}%'.format(date1, pctWater_roi_1*100)),
#         ('{0}: {1:2.1f}%'.format(date2, pctWater_roi_2*100))))
        
#     ax.text(0.78, 0.1, textstr, transform=ax.transAxes, fontsize=12,
#             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', edgecolor='none', alpha=0.8))



# Function to plot the difference in water for modified Devries (raster_1) and GFM (raster_2) and UNOSAT (raster_3), performs a raster sum to identify agreement / disagreement in water extent
def map_water_difference_vsGFMvsUNOSAT(raster_1, raster_2, raster_3, date, orig_crs, labels_crs, rounding_x, rounding_y, ax, cmap, cols_to_use, norm, raster_mask, districts, upazillas):

    # For this version, we just care about agreement / disagreement, not which one agrees with what - so just three colors = 1, 2 or 3 models have water
#     # The sum values we want therefore are:
#     # Single model = 1 or 2 or 4
#     # Two models = 3 or 5 or 6
#     # Three models = 7
    
#     # Set the water values in Devries to 1, and in GFM to 2
#     water_Devries = copy.copy(raster_1)
#     water_Devries.values[water_Devries.values<0]=0
#     water_Devries.values[water_Devries.values==1]=1

#     water_GFM = copy.copy(raster_2)
#     water_GFM.values[water_GFM.values<0]=0
#     water_GFM.values[water_GFM.values==1]=2
    
#     water_UNOSAT = copy.copy(raster_3)
#     water_UNOSAT.values[water_UNOSAT.values<0]=0
#     water_UNOSAT.values[water_UNOSAT.values==1]=4
    
#     # Sum the water values from the three rasters, 1/2/4, 1=Dev, 2=GFM, 3=Dev+GFM, 4=UNOSAT, 5=Dev+UNOSAT, 6=GFM+UNOSAT, 7=Dev+GFM+UNOSAT(all three)
#     water_sum = np.add(water_Devries, water_GFM, water_UNOSAT)

#     # Apply the mask to the water_sum raster by performing an array sum operation
#     water_sum_roi = np.add(water_sum, raster_mask)

#     # Get the number of pixels of water for each date for the ROI only (pixel values in ROI have been shifted by +10 from original values)
#     water_sum_roi_vals = list(chain.from_iterable(list(chain.from_iterable(water_sum_roi.values.tolist()))))
#     counter_roi = collections.Counter(water_sum_roi_vals)
#     total_roi = (counter_roi[10]+counter_roi[11]+counter_roi[12]+counter_roi[13]+counter_roi[14]+counter_roi[15]+counter_roi[16]+counter_roi[17])
#     pctWater_1model = (counter_roi[11]+counter_roi[12]+counter_roi[14]) / total_roi
#     pctWater_2model = (counter_roi[13]+counter_roi[15]+counter_roi[16]) / total_roi
#     pctWater_3model = (counter_roi[17]) / total_roi

# DO A SIMPLE VERSION INSTEAD, EASIER FOR PLOTTING WITH COLORMAP
    
    # Set the water values in each model to 1
    water_Devries = copy.deepcopy(raster_1)
    water_Devries.values[water_Devries.values<0]=0
    water_Devries.values[water_Devries.values==1]=1

    water_GFM = copy.deepcopy(raster_2)
    water_GFM.values[water_GFM.values<0]=0
    water_GFM.values[water_GFM.values==1]=1
    
    water_UNOSAT = copy.deepcopy(raster_3)
    water_UNOSAT.values[water_UNOSAT.values<0]=0
    water_UNOSAT.values[water_UNOSAT.values==1]=1
    
    # Sum the water values from the three rasters, 1=1model, 2=2models, 3=3models
    water_sum = np.add(water_Devries, water_GFM)
    water_sum = np.add(water_sum, water_UNOSAT)

#     # Apply the mask to the water_sum raster by performing an array sum operation
#     water_sum_roi = np.add(water_sum, raster_mask)

#     # Get the number of pixels of water for each date for the ROI only (pixel values in ROI have been shifted by +10 from original values)
#     water_sum_roi_vals = list(chain.from_iterable(list(chain.from_iterable(water_sum_roi.values.tolist()))))
#     counter_roi = collections.Counter(water_sum_roi_vals)
#     total_roi = (counter_roi[10]+counter_roi[11]+counter_roi[12]+counter_roi[13])
#     pctWater_1model = counter_roi[11] / total_roi
#     pctWater_2model = counter_roi[12] / total_roi
#     pctWater_3model = counter_roi[13] / total_roi   
    
#     print(counter_roi)
#     print('ROI % water for 1 model only: {0:2.1f}%'.format(pctWater_1model*100))
#     print('ROI % water for 2 models: {0:2.1f}%'.format(pctWater_2model*100))
#     print('ROI % water for 3 models: {0:2.1f}%'.format(pctWater_3model*100))

    # Create plot
    bbox = prep_raster.get_bbox_rioxarray(water_sum)
    raster_plot = ax.imshow(np.squeeze(water_sum), cmap=cmap, norm=norm, 
                            extent = [bbox.left, bbox.right, bbox.bottom, bbox.top], interpolation='none')

    # Add the Sylhet District outline
    #sylhet_dist.to_crs(epsg=32646).boundary.plot(ax=ax, color = 'yellow', linewidth=2, label='Sylhet District', zorder=3)
    districts.to_crs(epsg=32646).boundary.plot(ax=ax, color = 'yellow', linewidth=2, zorder=2)
    #sylhet_upa.to_crs(epsg=32646).plot(ax=ax, color = 'yellow', alpha=0.3, linewidth=1, label='Sylhet Upazilla', zorder=3)
    upazillas.to_crs(epsg=32646).boundary.plot(ax=ax, color = 'yellow', linewidth=1, zorder=2)

    # # Set title and legend
    ax.set_title(str(date))
    # legend_labels = {cols_to_use[1]: 'Background',
    #                  cols_to_use[2]: '1 model',
    #                  cols_to_use[3]: '2 models',
    #                  cols_to_use[4]: '3 models'}
    # patches = [Patch(color=color, label=label) for color, label in legend_labels.items()]
    # ax.legend(handles=patches, #bbox_to_anchor=(1.35, 1),
    #           facecolor='white')

    # Set the ticks and labels
    # Run the function to get the ticks and ticklabels
    # xticks_origCRS, yticks_origCRS, xticks_latlon, yticks_latlon = create_coord_ticks(bbox, orig_crs, labels_crs, rounding_x, rounding_y)
    # ax.set_xticks(xticks_origCRS, ['{:.1f}'.format(x) +'\N{DEGREE SIGN}E' for x in yticks_latlon])
    # ax.set_yticks(yticks_origCRS, ['{:.1f}'.format(x) +'\N{DEGREE SIGN}N' for x in xticks_latlon])
    ax.set_xlim(bbox.left, bbox.right)
    ax.set_ylim(bbox.bottom, bbox.top)
    ax.tick_params(top=False, bottom=False, left=False, right=False,
            labelleft=False, labelbottom=False)

    
    
    
# Function to plot the water across two dates, performs a raster sum to identify agreement / disagreement in water extent
def map_water_agree_subplot(source, raster_1, raster_2, name1, name2, orig_crs, labels_crs, rounding_x, rounding_y, ax, cmap, cols_to_use, norm, raster_mask, districts, upazillas):
    
    # Set the water values in source1 to 1, and in name2 to 2
    water_1 = copy.deepcopy(raster_1)
    water_1.values[water_1.values<0]=255
    water_1.values[water_1.values>1]=255
    water_1.values[water_1.values==1]=1

    water_2 = copy.deepcopy(raster_2)
    water_2.values[water_2.values<0]=255
    water_2.values[water_2.values>1]=255
    water_2.values[water_2.values==1]=1
    
    # Sum the water values from the two rasters, 0=background, 1=water one image, 2=water two images, set all else to nodataval
    water_sum = np.add(water_1, water_2)

    # Apply the mask to the water_sum raster by performing an array sum operation
    water_sum_roi = np.add(water_sum, raster_mask)

    # Get the number of pixels of water for each date for the ROI only (pixel values in ROI have been shifted by +10 from original values)
    water_sum_roi_vals = list(chain.from_iterable(list(chain.from_iterable(water_sum_roi.values.tolist()))))
    counter_roi = collections.Counter(water_sum_roi_vals) # 10 = none, 11 = water one image, 12 = water two images
    print(counter_roi)
    pctWater_agree_water = (counter_roi[12]) / (counter_roi[11] + counter_roi[12])
    pctWater_disagree_water = (counter_roi[11]) / (counter_roi[11] + counter_roi[12])
    pctWater_agree_all = (counter_roi[10] + counter_roi[12]) / (counter_roi[10] + counter_roi[11] + counter_roi[12])

    # Create plot
    bbox = prep_raster.get_bbox_rioxarray(water_sum)
    raster_plot = ax.imshow(np.squeeze(water_sum), cmap=cmap, norm=norm, 
                            extent = [bbox.left, bbox.right, bbox.bottom, bbox.top], interpolation='none', zorder=1)

    # Add the Sylhet District outline
    districts.to_crs(epsg=32646).boundary.plot(ax=ax, color = 'yellow', linewidth=1, zorder=2)
    upazillas.to_crs(epsg=32646).boundary.plot(ax=ax, color = 'yellow', linewidth=0.5, zorder=2)

    # Set title and legend
    ax.set_title(str(source), fontsize=20)

    # Set the ticks and labels
    # Run the function to get the ticks and ticklabels
    # xticks_origCRS, yticks_origCRS, xticks_latlon, yticks_latlon = create_coord_ticks(bbox, orig_crs, labels_crs, rounding_x, rounding_y)
    # ax.set_xticks(xticks_origCRS, ['{:.1f}'.format(x) +'\N{DEGREE SIGN}E' for x in yticks_latlon])
    # ax.set_yticks(yticks_origCRS, ['{:.1f}'.format(x) +'\N{DEGREE SIGN}N' for x in xticks_latlon])
    ax.set_xlim(bbox.left, bbox.right)
    ax.set_ylim(bbox.bottom, bbox.top)
    ax.tick_params(top=False, bottom=False, left=False, right=False,
                labelleft=False, labelbottom=False)
    
    # Add text box with the fractional area inundated
    textstr = '\n'.join((('Agreement water: {0:2.1f}%'.format(pctWater_agree_water*100)),
        ('Agreement all: {0:2.1f}%'.format(pctWater_agree_all*100))))
        
    ax.text(0.62, 0.12, textstr, transform=ax.transAxes, fontsize=16,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)) #edgecolor='none'

    
# Function to plot the water across two dates, performs a raster sum to identify agreement / disagreement in water extent
def map_waterMatchedDiff_subplot(source, raster_1, raster_2, bboxRaster, name1, name2, orig_crs, labels_crs, 
                                 rounding_x, rounding_y, ax, cmap, vmin, vmax, raster_mask, districts, upazillas):

    # Set the water values in source1 to 1, and in name2 to 2
    water_1 = copy.deepcopy(raster_1)
    water_2 = copy.deepcopy(raster_2)

    # Subtract the water values from the two rasters to get the absolute FIA difference
    water_diff = np.subtract(water_1, water_2) #subtracts water_2 from water_1, values could range from -1 to 1
    water_diff.values[np.logical_or(water_diff.values<-1,water_diff.values>1)]=None

    # Apply the mask to the water_sum raster by performing an array sum operation
    water_diff_roi = np.add(water_diff, raster_mask)

    # Get the mean FIA difference (only take values in the ROI i.e. values between 9 and 11
    water_diff_roi_mean = np.mean(water_diff_roi.values[np.logical_and(water_diff_roi.values>=9, water_diff_roi.values<=11)]) -10
    water_diff_roi_medi = np.median(water_diff_roi.values[np.logical_and(water_diff_roi.values>=9, water_diff_roi.values<=11)]) -10
    
    # Create plot
    bbox = prep_raster.get_bbox_rioxarray(water_diff)
    raster_plot = ax.imshow(np.squeeze(water_diff), cmap=cmap, vmin=vmin, vmax=vmax,
                        extent = [bbox.left, bbox.right, bbox.bottom, bbox.top], interpolation='none', zorder=1)

    # Add the Sylhet District outline
    districts.to_crs(epsg=32646).boundary.plot(ax=ax, color = 'yellow', linewidth=1, zorder=2)
    upazillas.to_crs(epsg=32646).boundary.plot(ax=ax, color = 'yellow', linewidth=0.5, zorder=2)

    # Set title and legend
    ax.set_title(str(source), fontsize=20)

    # Set the ticks and labels
    # Run the function to get the ticks and ticklabels
    # xticks_origCRS, yticks_origCRS, xticks_latlon, yticks_latlon = create_coord_ticks(bbox, orig_crs, labels_crs, rounding_x, rounding_y)
    # ax.set_xticks(xticks_origCRS, ['{:.1f}'.format(x) +'\N{DEGREE SIGN}E' for x in yticks_latlon])
    # ax.set_yticks(yticks_origCRS, ['{:.1f}'.format(x) +'\N{DEGREE SIGN}N' for x in xticks_latlon])
    ax.set_xlim(bbox.left, bbox.right)
    ax.set_ylim(bbox.bottom, bbox.top)
    ax.tick_params(top=False, bottom=False, left=False, right=False,
                labelleft=False, labelbottom=False)

    # Add text box with the fractional area inundated
    textstr = '\n'.join((('Mean FIA difference: {0:0.2f}'.format(water_diff_roi_mean)),
        ('Median FIA difference: {0:0.2f}'.format(water_diff_roi_medi))))

    ax.text(0.55, 0.12, textstr, transform=ax.transAxes, fontsize=16,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)) #edgecolor='none'
    
    # Add colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.2)
    cbar = plt.colorbar(raster_plot, cax=cax, ticks=[-1, 0, 1])
    ticklabs = cbar.ax.get_yticklabels()
    cbar.ax.set_yticklabels(ticklabs, fontsize=16)
    
    
# Function to get the map data for a given image and matched to a given label
def get_mapToPlot_wLabel(labelFile, mapRasterFile):
    
    # Get the label date and tileID
    labelDate=labelFile.stem.split('_')[0]
    labelTile='_'.join(labelFile.stem.split('_')[1:3])
    
    # Open the label raster
    labelRaster = prep_raster.load_label_for_comparison(str(labelFile))
    
    # Open the map to compare, and reproject and match the extent
    mapRaster = prep_raster.load_raster_for_comparison(str(mapRasterFile), 0)
    mapRasterReproj = prep_raster.reproj_match_raster(mapRaster, labelRaster)
       
    return labelRaster, mapRasterReproj
    
# Function to plot 5 by 1 subplot showing a single tile with label and the maps and their IOU scores    
def plotCompare_wLabel_wFFWC(labelSelect, labelFiles, 
                       compareDate, Devries_mosaics, GFM_mosaics, IMPACT_mosaics,
                       FFWCthresh, FFWCFiles,
                       accuracyMaps, accuracyFFWC,
                       cmap_to_use, cols_to_use, norm):
    
    # Select the label and get its date, get the label file
    labelDate=labelSelect.split('_')[0]
    labelDateFile = [file for file in labelFiles if labelSelect in str(file)][0]
    
    # Get the Devries, GFM and NASAIMPACT rasters to comapre, from the closest date
    DevriesCompare = [file for file in Devries_mosaics if compareDate in str(file)][0]
    GFMCompare = [file for file in GFM_mosaics if str(datetime.strptime(compareDate, '%Y%m%d').date()).replace('-','_') in str(file)][0]
    IMPACTCompare = [file for file in IMPACT_mosaics if compareDate in str(file)][0]
    
    # Get the FFWC map to compare
    FFWCCompare = [file for file in FFWCFiles if str(FFWCthresh) in str(file)][0]
    FFWCdate = FFWCCompare.stem.split('_')[0]
    
    # Get the rasters in the correct format for plotting
    labelPlot, DevriesPlot = get_mapToPlot_wLabel(labelDateFile, DevriesCompare)
    labelPlot, GFMPlot = get_mapToPlot_wLabel(labelDateFile, GFMCompare)
    labelPlot, IMPACTPlot = get_mapToPlot_wLabel(labelDateFile, IMPACTCompare)
    labelPlot, FFWCPlot = get_mapToPlot_wLabel(labelDateFile, FFWCCompare)
    
    # Get the IOU values to add as plot annotation
    mapsIOU = [accuracyMaps[accuracyMaps.label==labelSelect]['IOU'].values[0]]+[accuracyMaps[accuracyMaps.label==labelSelect]['IOU'].values[2]]+[accuracyMaps[accuracyMaps.label==labelSelect]['IOU'].values[1]]
    FFWCIOU = accuracyFFWC[np.logical_and(accuracyFFWC.label==labelSelect, accuracyFFWC.DepthThreshold==float(FFWCthresh))]['IOU'].values[0]
    rasterIOU = [''] + mapsIOU + [FFWCIOU]
    
    # Organize the rasters for plotting
    rasterDataPlot = [labelPlot, DevriesPlot, IMPACTPlot, GFMPlot, FFWCPlot]
    rasterDates = [labelDate, compareDate, compareDate, compareDate, FFWCdate]
    rasterDates = [str(datetime.strptime(item, '%Y%m%d').date()) for item in rasterDates]
    rasterNames = ['Label', 'Thomas et al.', 'Paul & Ganju', 'GFM', 'FFWC']
    rasterThresh = ['','','','',''] #str(FFWCthresh)
    rasterStringsDF = pd.DataFrame(data=[rasterNames, rasterDates, rasterThresh, rasterIOU], index=['name','date','thresh', 'iou']).T
    rasterStringsDF['string']=rasterStringsDF['name'].map(str)+' '+rasterStringsDF['date'].map(str)+' '+rasterStringsDF['thresh'].map(str)+' \n IOU: '+rasterStringsDF['iou'].astype(str).str[0:5]
    rasterStringsDF['string'][0]=rasterStringsDF['string'][0].replace('IOU:','')
    
    fig, ax = plt.subplots(1, 5, figsize=(25,5))

    for i, axi in enumerate(ax.ravel()):

        axi.imshow(rasterDataPlot[i][0], norm=norm, cmap=cmap_to_use, interpolation='none') #vmin=0, vmax=1, 
        axi.set_title(rasterStringsDF.string.iloc[i], fontsize=18)
        axi.tick_params(top=False, bottom=False, left=False, right=False, labelleft=False, labelbottom=False)


    legend_labels =  {cols_to_use[0]: 'Background',
                     cols_to_use[1]: 'Water'}
    patches = [Patch(color=color, label=label) for color, label in legend_labels.items()]
    fig.legend(handles=patches, bbox_to_anchor=(0.5, 0),
                      loc='lower center', ncols=len(legend_labels),
                      facecolor='white', fontsize=18) 

    # fig.tight_layout()
    plt.show()
    
    
# Function to plot 4 by 1 subplot showing a single tile with label and the maps and their IOU scores    
def plotCompare_wLabel(labelSelect, labelFiles, 
                       compareDate, Devries_mosaics, GFM_mosaics, IMPACT_mosaics,
                       accuracyMaps, 
                       cmap_to_use, cols_to_use, norm):
    
    # Select the label and get its date, get the label file
    labelDate=labelSelect.split('_')[0]
    labelDateFile = [file for file in labelFiles if labelSelect in str(file)][0]
    
    # Get the Devries, GFM and NASAIMPACT rasters to comapre, from the closest date
    DevriesCompare = [file for file in Devries_mosaics if compareDate in str(file)][0]
    GFMCompare = [file for file in GFM_mosaics if str(datetime.strptime(compareDate, '%Y%m%d').date()).replace('-','_') in str(file)][0]
    IMPACTCompare = [file for file in IMPACT_mosaics if compareDate in str(file)][0]
    
    
    # Get the rasters in the correct format for plotting
    labelPlot, DevriesPlot = get_mapToPlot_wLabel(labelDateFile, DevriesCompare)
    labelPlot, GFMPlot = get_mapToPlot_wLabel(labelDateFile, GFMCompare)
    labelPlot, IMPACTPlot = get_mapToPlot_wLabel(labelDateFile, IMPACTCompare)
    
    # Get the IOU values to add as plot annotation
    mapsIOU = [accuracyMaps[accuracyMaps.label==labelSelect]['IOU'].values[0]]+[accuracyMaps[accuracyMaps.label==labelSelect]['IOU'].values[2]]+[accuracyMaps[accuracyMaps.label==labelSelect]['IOU'].values[1]]
    rasterIOU = [''] + mapsIOU
    
    # Organize the rasters for plotting
    rasterDataPlot = [labelPlot, DevriesPlot, IMPACTPlot, GFMPlot]
    rasterDates = [labelDate, compareDate, compareDate, compareDate]
    rasterDates = [str(datetime.strptime(item, '%Y%m%d').date()) for item in rasterDates]
    rasterNames = ['Label', 'Thomas et al.', 'Paul & Ganju', 'GFM']
    rasterThresh = ['','','',''] #str(FFWCthresh)
    rasterStringsDF = pd.DataFrame(data=[rasterNames, rasterDates, rasterIOU], index=['name','date','iou']).T
    rasterStringsDF['string']=rasterStringsDF['name'].map(str)+' '+rasterStringsDF['date'].map(str)+' \n IOU: '+rasterStringsDF['iou'].astype(str).str[0:5]
    rasterStringsDF['string'][0]=rasterStringsDF['string'][0].replace('IOU:','')
    
    fig, ax = plt.subplots(1, 4, figsize=(20,5))

    for i, axi in enumerate(ax.ravel()):

        axi.imshow(rasterDataPlot[i][0], norm=norm, cmap=cmap_to_use, interpolation='none') #vmin=0, vmax=1, 
        axi.set_title(rasterStringsDF.string.iloc[i], fontsize=18)
        axi.tick_params(top=False, bottom=False, left=False, right=False, labelleft=False, labelbottom=False)


    legend_labels =  {cols_to_use[0]: 'Background',
                     cols_to_use[1]: 'Water'}
    patches = [Patch(color=color, label=label) for color, label in legend_labels.items()]
    fig.legend(handles=patches, bbox_to_anchor=(0.5, 0),
                      loc='lower center', ncols=len(legend_labels),
                      facecolor='white', fontsize=18) 

    # fig.tight_layout()
    plt.show()
    
    
# Function to plot 4 by 1 subplot showing a single tile with label and the maps and their IOU scores   
# Plots in confusion matrix different colors for TP, FP, TN, FN etc.
def plotCompare_wLabel_confusMat(labelSelect, labelFiles, 
                       compareDate, Devries_mosaics, GFM_mosaics, IMPACT_mosaics,
                       accuracyMaps, 
                       cmap_to_use, cols_to_use, norm,
                       cmap_label, cols_label, norm_label):
    
    # Select the label and get its date, get the label file
    labelDate=labelSelect.split('_')[0]
    labelDateFile = [file for file in labelFiles if labelSelect in str(file)][0]
    
    # Get the Devries, GFM and NASAIMPACT rasters to comapre, from the closest date
    DevriesCompare = [file for file in Devries_mosaics if compareDate in str(file)][0]
    GFMCompare = [file for file in GFM_mosaics if str(datetime.strptime(compareDate, '%Y%m%d').date()).replace('-','_') in str(file)][0]
    IMPACTCompare = [file for file in IMPACT_mosaics if compareDate in str(file)][0]
    
    
    # Get the rasters in the correct format for plotting
    labelPlot, DevriesPlot = get_mapToPlot_wLabel(labelDateFile, DevriesCompare)
    labelPlot, GFMPlot = get_mapToPlot_wLabel(labelDateFile, GFMCompare)
    labelPlot, IMPACTPlot = get_mapToPlot_wLabel(labelDateFile, IMPACTCompare)
    
    # Get the IOU values to add as plot annotation
    mapsIOU = [accuracyMaps[accuracyMaps.label==labelSelect]['IOU'].values[0]]+[accuracyMaps[accuracyMaps.label==labelSelect]['IOU'].values[2]]+[accuracyMaps[accuracyMaps.label==labelSelect]['IOU'].values[1]]
    rasterIOU = [''] + mapsIOU
    
    # Organize the rasters for plotting
    rasterDataPlot = [labelPlot, DevriesPlot, IMPACTPlot, GFMPlot]
    rasterDates = [labelDate, compareDate, compareDate, compareDate]
    rasterDates = [str(datetime.strptime(item, '%Y%m%d').date()) for item in rasterDates]
    rasterNames = ['Label', 'Thomas et al.', 'Paul & Ganju', 'GFM']
    rasterThresh = ['','','',''] #str(FFWCthresh)
    rasterStringsDF = pd.DataFrame(data=[rasterNames, rasterDates, rasterIOU], index=['name','date','iou']).T
    rasterStringsDF['string']=rasterStringsDF['name'].map(str)+' '+rasterStringsDF['date'].map(str)+' \n IOU: '+rasterStringsDF['iou'].astype(str).str[0:5]
    rasterStringsDF['string'][0]=rasterStringsDF['string'][0].replace('IOU:','')
    
    fig, ax = plt.subplots(1, 4, figsize=(14,4), layout='constrained') #figsize=(20,5)

    for i, axi in enumerate(ax.ravel()):
        
        
        # Label
        if i==0:
            axi.imshow(rasterDataPlot[i][0], norm=norm_label, cmap=cmap_label, interpolation='none') #vmin=0, vmax=1, 
            axi.set_title(rasterStringsDF.string.iloc[i], fontsize=18)
            axi.tick_params(top=False, bottom=False, left=False, right=False, labelleft=False, labelbottom=False)
        
        # Others
        else:
            # Get the raster and compare with the label, get TP, TN, FP, FNs
            raster=rasterDataPlot[i][0]          
            label=rasterDataPlot[0][0]
            
            raster_CM=copy.deepcopy(raster)
            raster_CM.values[np.logical_and(raster.values==1,label.values==1)]=0 # TP
            raster_CM.values[np.logical_and(raster.values==0,label.values==0)]=1 # TN
            raster_CM.values[np.logical_and(raster.values==1,label.values==0)]=2 # FP
            raster_CM.values[np.logical_and(raster.values==0,label.values==1)]=3 # FN

            axi.imshow(raster_CM, norm=norm, cmap=cmap_to_use, interpolation='none') #vmin=0, vmax=1, 
            axi.set_title(rasterStringsDF.string.iloc[i], fontsize=18)
            axi.tick_params(top=False, bottom=False, left=False, right=False, labelleft=False, labelbottom=False)


    legend_labels =  {cols_label[0]: 'Background',
                     cols_label[1]: 'Water',
                     cols_to_use[0]: 'True Positive',
                     cols_to_use[1]: 'True Negative',
                     cols_to_use[2]: 'False Positive',
                     cols_to_use[3]: 'False Negative'}
    patches = [Patch(color=color, label=label) for color, label in legend_labels.items()]
    fig.legend(handles=patches, bbox_to_anchor=(0.5, -0.08),
                      loc='center', ncols=len(legend_labels),
                      facecolor='white', fontsize=16) 

    # fig.tight_layout()
    plt.show()

    
    
    
    
    
# Function to plot a scatter and regression of FIA values from two rasters matched to Fusion, as a subplot    

def subplotFIAScatter(ax, raster1, raster2):

    # Define x and y for plottin and regression
    x = raster1[0].values.ravel()
    y = raster2[0].values.ravel()
    x[x==255]=np.nan
    y[y==255]=np.nan

    # Calculate correlation, slope, int stats and record
    mask = ~np.isnan(x) & ~np.isnan(y)
    stats_vals = scipy.stats.linregress(x[mask], y[mask]) #slope, intercept, r_value, p_value, std_err

    # Plot scatter points
    ax.scatter(x, y, marker='o', s=1, color = 'b', alpha=0.3)

    # Plot linear trend
    ax.plot(x, stats_vals[0]*x+stats_vals[1], color = 'r', linewidth=2, zorder=3)
    ax.plot([0,1], [0,1], color = 'k', linewidth=2, linestyle='--', zorder=2)

    # Add annotation with slope and correlation values
    slope = stats_vals[0]
    correl = stats_vals[2]**2
    stats_string = 'm: {0:1.2f}, $R^2$: {1:1.2f}'.format(slope, correl)

    # Set the yaxis, gridlines and limits
    ax.set_xticks([0.05,0.95], ['0','1'])
    ax.set_yticks([0.05,0.95], ['','1'])
    ax.tick_params(axis='both', which='major', direction='in', pad=-20, length=0, labelsize=20, labelcolor='w')
    ax.set_ylim(0,1)
    ax.set_xlim(0,1)
    
    xtick_labels = ax.get_xticklabels()
    ytick_labels = ax.get_yticklabels()
    white_effect = pe.Stroke(linewidth=3, foreground='k')
    for label in xtick_labels:
        label.set_path_effects([white_effect,pe.Normal()])
    for label in ytick_labels:
        label.set_path_effects([white_effect,pe.Normal()])

    # Annotate the slope and correlation values
    ax.text(0.05, 0.70, stats_string, transform=ax.transAxes, fontsize=20, color='k',
            verticalalignment='bottom', bbox=dict(boxstyle='round', facecolor='white', edgecolor='lightgray', alpha=0.8))
    
    
# Function to plot the water across two dates, performs a raster sum to identify agreement / disagreement in water extent - for the pairs plot
def subplotFIADiffMap(raster_1, raster_2, ax, cmap, vmin, vmax, districts, upazillas):

    # Set the water values in source1 to 1, and in name2 to 2
    water_1 = copy.deepcopy(raster_1[0])
    water_2 = copy.deepcopy(raster_2[0])

    # Subtract the water values from the two rasters to get the absolute FIA difference
    water_diff = np.subtract(water_1, water_2) #subtracts water_2 from water_1, values could range from -1 to 1
    water_diff.values[np.logical_or(water_diff.values<-1,water_diff.values>1)]=None

#     # Apply the mask to the water_sum raster by performing an array sum operation
#     water_diff_roi = np.add(water_diff, raster_mask)

#     # Get the mean FIA difference (only take values in the ROI i.e. values between 9 and 11
#     water_diff_roi_mean = np.mean(water_diff_roi.values[np.logical_and(water_diff_roi.values>=9, water_diff_roi.values<=11)]) -10
#     water_diff_roi_medi = np.median(water_diff_roi.values[np.logical_and(water_diff_roi.values>=9, water_diff_roi.values<=11)]) -10
    
    # Create plot
    bbox = prep_raster.get_bbox_rioxarray(water_diff)
    aspectRatio=water_diff.shape[1]/water_diff.shape[0]
    raster_plot = ax.imshow(np.squeeze(water_diff), cmap=cmap, vmin=vmin, vmax=vmax, aspect=aspectRatio,
                        extent = [bbox.left, bbox.right, bbox.bottom, bbox.top], interpolation='none', zorder=1)

    # Add the Sylhet District outline
    districts.to_crs(epsg=32646).boundary.plot(ax=ax, color = 'k', linewidth=2, zorder=2)
    ax.set_aspect(aspectRatio)
    # upazillas.to_crs(epsg=32646).boundary.plot(ax=ax, color = 'k', linewidth=0.5, zorder=2)

    # Set the ticks and labels
    ax.set_xlim(bbox.left, bbox.right)
    ax.set_ylim(bbox.bottom, bbox.top)
    ax.tick_params(top=False, bottom=False, left=False, right=False,
                labelleft=False, labelbottom=False)

#     # Add text box with the fractional area inundated
#     textstr = '\n'.join((('Mean FIA difference: {0:0.2f}'.format(water_diff_roi_mean)),
#         ('Median FIA difference: {0:0.2f}'.format(water_diff_roi_medi))))

#     ax.text(0.55, 0.12, textstr, transform=ax.transAxes, fontsize=16,
#             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)) #edgecolor='none'
    
    # Return the imshow plot to use for adding a colorbar
    return raster_plot


# Function to plot the water one single date, - for the pairs plot
def subplotFIAMap(raster_1, orig_crs, labels_crs, 
                                 rounding_x, rounding_y, ax, cmap, vmin, vmax, districts, upazillas, ticks=False, sbar=False):

    # Set the water values in source1 to 1, and in name2 to 2
    water_1 = copy.deepcopy(raster_1[0])

    # Subtract the water values from the two rasters to get the absolute FIA difference
    water_1.values[np.logical_or(water_1.values<0,water_1.values>1)]=None

    # Create plot
    bbox = prep_raster.get_bbox_rioxarray(water_1)
    aspectRatio=water_1.shape[1]/water_1.shape[0]
    raster_plot = ax.imshow(np.squeeze(water_1), cmap=cmap, vmin=vmin, vmax=vmax, aspect=aspectRatio,
                        extent = [bbox.left, bbox.right, bbox.bottom, bbox.top], interpolation='none', zorder=1)

    # Add the Sylhet District outline
    districts.to_crs(epsg=32646).boundary.plot(ax=ax, color = 'k', linewidth=2, zorder=1)
    ax.set_aspect(aspectRatio)
    # upazillas.to_crs(epsg=32646).boundary.plot(ax=ax, color = 'k', linewidth=0.5, zorder=1)

    # Set the ticks and labels
    # Run the function to get the ticks and ticklabels
    ax.set_xlim(bbox.left, bbox.right)
    ax.set_ylim(bbox.bottom, bbox.top)
    if ticks==True:
        xticks_origCRS, yticks_origCRS, xticks_latlon, yticks_latlon = create_coord_ticks(bbox, orig_crs, labels_crs, rounding_x, rounding_y)
        xkeep = [1,3]
        ykeep = [0,2]
        xticks_origCRS, yticks_origCRS, xticks_latlon, yticks_latlon = [xticks_origCRS[i] for i in xkeep], [yticks_origCRS[i] for i in ykeep], [xticks_latlon[i] for i in xkeep], [yticks_latlon[i] for i in ykeep]
        ax.set_xticks(xticks_origCRS, ['{:.1f}'.format(x) +'\N{DEGREE SIGN}E' for x in yticks_latlon])
        ax.set_yticks(yticks_origCRS, ['{:.1f}'.format(x) +'\N{DEGREE SIGN}N' for x in xticks_latlon])
        ax.tick_params(axis='x', which='major', direction='in', pad=-35, width=3, length=10, color='k',
                       labelsize=20, labelcolor='w', grid_color='k', grid_linestyle='--', grid_alpha=0.8, zorder=3)
        ax.tick_params(axis='y', which='major', direction='in', pad=-75, width=3, length=10, color='k',
                       labelsize=20, labelcolor='w', grid_color='k', grid_linestyle='--', grid_alpha=0.8, zorder=3)

        # Get the xtick labels
        xtick_labels = ax.get_xticklabels()
        ytick_labels = ax.get_yticklabels()
        white_effect = pe.Stroke(linewidth=3, foreground='k')
        for label in xtick_labels:
            label.set_path_effects([white_effect,pe.Normal()])
        for label in ytick_labels:
            label.set_path_effects([white_effect,pe.Normal()])

    else:
        ax.tick_params(top=False, bottom=False, left=False, right=False,
                labelleft=False, labelbottom=False)
        
    # Add a scalebar
    if sbar==True:
        fontprops = fm.FontProperties(size=20)
        scalebar = AnchoredSizeBar(ax.transData,
                               20000, '20 km', 'upper right', 
                               pad=0.7, borderpad=0.05,
                               color='black', frameon=False,
                               size_vertical=20000/20,
                               fontproperties=fontprops)
        ax.add_artist(scalebar)

   
    # Return the imshow plot to use for adding a colorbar
    return raster_plot