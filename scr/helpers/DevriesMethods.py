import ee

# elevation = ee.Image("CGIAR/SRTM90_V4")


def calc_basemean (x, start, end, mode, direction):
    return x \
    .filter(ee.Filter.equals("orbitProperties_pass", direction)) \
    .filter(ee.Filter.equals("instrumentMode", mode)) \
    .filterDate(start, end) \
    .mean()


# baseline (sd)
def calc_basesd (x, start, end, mode, direction):
    return x \
    .filter(ee.Filter.equals("orbitProperties_pass", direction)) \
    .filter(ee.Filter.equals("instrumentMode", mode)) \
    .filterDate(start, end) \
    .reduce(ee.Reducer.stdDev()) \
    .rename(["VV", "VH", "angle"])


# anomaly
def calc_anomaly(x, start, end, mode, direction):
    basemean = calc_basemean(x, start, end, mode, direction)

    return x\
    .filter(ee.Filter.equals('orbitProperties_pass', direction))\
    .filter(ee.Filter.equals('instrumentMode', mode))\
    .map(lambda y:y\
        .subtract(basemean)\
        .set({'system:time_start': y.get('system:time_start'),\
        'orbitProperties_pass': direction,\
        'orbitNumber_start': y.get('orbitNumber_start')})\
        )

def calc_zscore(x, start, end, mode, direction):
    anom = calc_anomaly(x, start, end, mode, direction)
    basesd = calc_basesd(x, start, end, mode, direction)
    return anom.map(lambda xElement: xElement.divide(basesd).set('system:time_start', xElement.get('system:time_start')).set('orbitNumber_start', xElement.get('orbitNumber_start')).set('orbitProperties_pass', direction))

def add_dt_str (image):
    date = ee.Date(image.get("system:time_start"))
    return image.set("Date", date.format()).copyProperties(image)

def addIndexBands(image):

    # band5 ndwi
    ndwi2 = image.normalizedDifference(['sur_refl_b04','sur_refl_b05']).rename(['NDWI2'])

    # lswi
    lswi = image.normalizedDifference(['sur_refl_b02','sur_refl_b05']).rename(['LSWI'])

    # b1/b2
    b1_div_b2 = image.select('sur_refl_b01').divide(image.select('sur_refl_b02')).rename(['B1/B2'])

    # ndvi
    ndvi = image.normalizedDifference(['sur_refl_b02','sur_refl_b01']).rename(['NDVI'])

    elevation = ee.ImageCollection("projects/sat-io/open-datasets/FABDEM").mosaic()
  
    return image.addBands(ndwi2).addBands(lswi).addBands(b1_div_b2).addBands(ndvi).addBands(elevation.rename('elevation'))

