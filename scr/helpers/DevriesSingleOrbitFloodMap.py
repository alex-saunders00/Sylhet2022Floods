import Devries.DevriesMethods as dv
import ee


def makeBoolFloodMap (geometry, basestart,  baseend, orbitNumber,  kernelRadius,  kernelThreshold,  addThresholding,  vhThreshold,  z_th,  smooth ):
    
    filters = [
        ee.Filter.listContains("transmitterReceiverPolarisation", "VV"),\
        ee.Filter.listContains("transmitterReceiverPolarisation", "VH"),\
        ee.Filter.equals("instrumentMode", "IW"),\
        ee.Filter.bounds(geometry),
    ]

    # Load S1 collection with filters
    s1 = ee.ImageCollection("COPERNICUS/S1_GRD").filter(filters).filter(ee.Filter.eq('orbitNumber_start', int(orbitNumber)))
    s1ProjImage = s1.first().select('VH')
    s1_combined = ee.ImageCollection("COPERNICUS/S1_GRD").filter(filters).filterDate(basestart,baseend)\
        .merge(s1)\
        .map(lambda image: image.updateMask(image.select('VV').gt(-45)))

    s1_asc = s1_combined.filter(ee.Filter.equals('orbitProperties_pass', 'ASCENDING'));
    s1_dsc = s1_combined.filter(ee.Filter.equals('orbitProperties_pass', 'DESCENDING'));

    z_iwasc = dv.calc_zscore(s1_asc, basestart, baseend, 'IW', 'ASCENDING');
    z_iwdsc = dv.calc_zscore(s1_dsc, basestart, baseend, 'IW', 'DESCENDING');
    z = ee.ImageCollection(z_iwasc.merge(z_iwdsc)).sort('system:time_start').map(dv.add_dt_str);
    
    z_study = z.filter(ee.Filter.eq('orbitNumber_start', int(orbitNumber)))

    #use medium con fidence and high (so vv OR vh)
    vvflag = z_study.mosaic().select('VV').lt(z_th)
    vhflag = z_study.mosaic().select('VH').lt(z_th)

    boolImage = vvflag.Or(vhflag)

    # add raw VH threshold
    if(addThresholding):
        vhraw = s1_combined.select('VH').mosaic().lt(vhThreshold);
        boolImage = boolImage.add(vhraw).gte(1);

    # smooth image if true
    if(smooth):
        boolImage = boolImage.setDefaultProjection(s1ProjImage.projection(), None, 10).focal_mean(radius=kernelRadius, units='meters')\
          .gt(kernelThreshold).setDefaultProjection(s1ProjImage.projection());

    mask = s1.mosaic().select('VV').unmask(-9999999999999).gt(-40).selfMask().mask()

    boolImageMasked = boolImage.mask(mask).unmask(-9999999999999)
    
    boolImageReduced = boolImageMasked.select(0)\
                        .setDefaultProjection(crs='EPSG:32646', scale=10)\
                        .reduceResolution(
                            reducer=ee.Reducer.mean(),\
                            maxPixels=3000\
                        )\
                        .rename(['s1_frac_flood'])
    return boolImageReduced

def getGTImageSingleOrbit(chip, basestart, baseend, orbitNumber, s1ImageDate):
    return makeBoolFloodMap(chip, basestart, baseend,\
                        orbitNumber,\
                        30,0.6,True,-22,-2, True)\
                        .set("s1date", s1ImageDate)
    