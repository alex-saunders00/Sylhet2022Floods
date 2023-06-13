import Devries.DevriesMethods as dv
import importlib
importlib.reload(dv)
import ee


def makeBoolFloodMap (geometry, basestart,  baseend,  targstart, targend,  kernelRadius,  kernelThreshold,  addThresholding,  vhThreshold,  z_th,  smooth ):
    targbefore = ee.Date(targstart).advance(-1, "month")

    # Define collection filters
    filters = [
        ee.Filter.listContains("transmitterReceiverPolarisation", "VV"),\
        ee.Filter.listContains("transmitterReceiverPolarisation", "VH"),\
        ee.Filter.equals("instrumentMode", "IW"),\
        ee.Filter.bounds(geometry),
    ]

    # Load S1 collection with filters
    s1 = ee.ImageCollection("COPERNICUS/S1_GRD").filter(filters)
    s1ProjImage = s1.filterDate(targstart,targend).first().select('VH')
    s1_combined = s1.filterDate(basestart,baseend)\
        .merge(s1.filterDate(targbefore,ee.Date(targend).advance(1,'day')))\
        .map(lambda image: image.updateMask(image.select('VV').gt(-45)))
    
    s1_asc = s1_combined.filter(ee.Filter.equals('orbitProperties_pass', 'ASCENDING'));
    s1_dsc = s1_combined.filter(ee.Filter.equals('orbitProperties_pass', 'DESCENDING'));

    z_iwasc = dv.calc_zscore(s1_asc, basestart, baseend, 'IW', 'ASCENDING');
    z_iwdsc = dv.calc_zscore(s1_dsc, basestart, baseend, 'IW', 'DESCENDING');
    z = ee.ImageCollection(z_iwasc.merge(z_iwdsc)).sort('system:time_start').map(dv.add_dt_str);
    
    z_study = z.filterDate(targstart, targend)

    #use medium confidence and high (so vv OR vh)
    vvflag = z_study.mosaic().select('VV').lt(z_th)
    vhflag = z_study.mosaic().select('VH').lt(z_th)

    boolImage = vvflag.Or(vhflag)
    # add raw VH threshold
    if(addThresholding):
        vhraw = s1_combined.select('VH').mosaic().lt(vhThreshold)
        boolImage = boolImage.add(vhraw).gte(1)

    # smooth image if true
    if(smooth):
        boolImage = boolImage.setDefaultProjection(s1ProjImage.projection(), None, 10).focal_mean(radius=kernelRadius, units='meters')\
          .gt(kernelThreshold).setDefaultProjection(s1ProjImage.projection())
    
    boolImage = boolImage.select(0)\
                        .setDefaultProjection(crs='EPSG:32646', scale=10)\
                        .reduceResolution(
                            reducer=ee.Reducer.mean(),\
                            maxPixels=3000\
                        )\
                        .rename(['s1_frac_flood'])
    return boolImage

def getModisGTCombined(modisCommand, chip, basestart, baseend, s1ImageDate):
    s1FloodImage = makeBoolFloodMap(chip, basestart, baseend,\
                        ee.Date(s1ImageDate).advance(-1,'month'), ee.Date(s1ImageDate).advance(1,'day'),\
                        30,0.6,True,-22,-2, True)
    
    modisImage = modisCommand.mosaic()
    modisImageWithS1Date = modisImage.set("s1date", s1ImageDate)
    modisImageWithS1Date = modisImageWithS1Date.addBands(modisImageWithS1Date.metadata("s1date"))

    combinedImage = dv.addIndexBands(modisImageWithS1Date).addBands(s1FloodImage)
    return combinedImage

def getGTImage(chip, basestart, baseend, s1ImageDate):
    return makeBoolFloodMap(chip, basestart, baseend,\
                        ee.Date(s1ImageDate).advance(-1,'month'), ee.Date(s1ImageDate).advance(1,'day'),\
                        30,0.6,True,-22,-2, True)\
                        .set("s1date", s1ImageDate)

def getModis(modisCommand, s1ImageDate):
    
    modisImage = modisCommand.mosaic()
    modisImageWithS1Date = modisImage.set("s1date", s1ImageDate)

    modisImageWithIndexes = dv.addIndexBands(modisImageWithS1Date)
    return modisImageWithIndexes