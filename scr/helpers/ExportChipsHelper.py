import pandas as pd
import ee
from pathlib import Path
from py_linq import Enumerable

def GetS1Dates(s1_raw, chip, globalStart, globalEnd):
    filters = [
        ee.Filter.listContains("transmitterReceiverPolarisation", "VV"),\
        ee.Filter.listContains("transmitterReceiverPolarisation", "VH"),\
        ee.Filter.equals("instrumentMode", "IW")
    ]

    S1Filtered = s1_raw.filterBounds(chip).filterDate(ee.Date(globalStart), ee.Date(globalEnd)).filter(filters)

    dates = S1Filtered.aggregate_array("system:time_start").getInfo()

    orbitNumber_start = S1Filtered.aggregate_array("orbitNumber_start").getInfo()
    sat = S1Filtered.aggregate_array("orbitProperties_pass").getInfo()
    platform = S1Filtered.aggregate_array("platform_number").getInfo()

    datesSorted, orbitNumber_startSorted, sat_sorted, platform_sorted = zip(*sorted(zip(dates, orbitNumber_start, sat, platform)))

    dfProperties = pd.DataFrame()
    dfProperties['orbitNumber'] = orbitNumber_startSorted
    dfProperties['sat'] = sat_sorted
    dfProperties['platform'] = platform_sorted
    dfProperties = dfProperties.drop_duplicates()

    df = pd.DataFrame()
    df['date'] = datesSorted
    df['orbitNumber'] = orbitNumber_startSorted

    dfGrouped = df.groupby('orbitNumber').agg({'date': ['min', 'max', 'median', 'count']})

    return dfGrouped, dfProperties

def createFineGrid(gridElement, scale=32*500, selectOnlyElementsInsideRoI=True):
    grid = gridElement.geometry().coveringGrid(proj=gridElement.geometry().projection(), scale=scale)
    
    if selectOnlyElementsInsideRoI:
        grid = grid.map(lambda subGridElement:subGridElement\
            .set('within', subGridElement.centroid(maxError=1).containedIn(gridElement.geometry())))\
            .filter(ee.Filter.eq('within', True))
    
    return grid

def getFileList(folder):
    Path.lsTif = lambda x: Enumerable(x.iterdir()).where(lambda p: p.suffix == '.tif').to_list()
    return [path.stem for path in folder.lsTif()]
