import pandas as pd
import ee
from pathlib import Path
from py_linq import Enumerable

import pandas as pd
import numpy as np
import ee
import geemap
import importlib
from pathlib import Path
import contextlib
from functools import partial

import multiprocessing

# import Devries.DevriesFloodMap as dvFM
# importlib.reload(dvFM)
# import Devries.DevriesMethods as dv
# importlib.reload(dv)
# import Devries.DevriesSingleOrbitFloodMap as dvSOFM
# importlib.reload(dvSOFM)
# import Devries.ExportChipsHelper as ExportHelper
# importlib.reload(ExportHelper)



# Some helper functions from the original authors
# https://github.com/sidgan/ETCI-2021-Competition-on-Flood-Detection/blob/main/notebooks/Ensemble_Inference.ipynb

def get_test_id(path):
    return path.split("_")[0] + "_" + path.split("_")[1]

def make_im_name(id, suffix):
    return id.split(".")[0] + f"_{suffix}.png"

def s1_to_rgb(vv_image, vh_image):
    ratio_image = np.clip(np.nan_to_num(vh_image/vv_image, 0), 0, 1)
    rgb_image = np.stack((vv_image, vh_image, 1-ratio_image), axis=2)
    return rgb_image


