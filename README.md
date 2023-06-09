# Sylhet2022Floods
This repository contains code used in the following paper, which was submitted to the Institute of Electrical and Electronics Engineers (IEEE) International Symposium on Geoscience and Remote Sensing (IGARSS) 2023 conference.

### A comparison of remote sensing approaches to assess the devastating May-June 2022 floods in Sylhet, Bangladesh

Saunders, A.<sup>1</sup>, Giezendanner, J.<sup>1</sup>, Tellman, B.<sup>1</sup>, Islam, A.<sup>1</sup>, Bhuyan, A.<sup>2</sup>, Islam, A.K.M.S.<sup>3</sup>

<sup>1</sup>Social Pixel Lab, School of Geography, Development and Environment, University of Arizona, USA
<sup>2</sup>Flood Forecasting and Warning Centre, Bangladesh Water Development Board
<sup>3</sup>Institute of Water and Flood Management, Bangladesh University of Engineering and Technology

### Description of the code
The code is divided into the following sections:
* **prep:** Data Preparation
* **analysis:** Data Analysis
* **results:** Results Analysis

The **prep** folder is split into a subsequent two folders for creating surface water maps using (1) the Thomas et al. (2023) ("local", non-machine learning) and (2) the Paul & Ganju (2021) ("global", machine learning) algorithms. Code for the Thomas et al. (2023) algorithm was adapted from the GitHub page https://github.com/mitchellthomas1/S1-Flood-Bangladesh. Code for the Paul & Ganju (2021) algorithm was adapted from the GitHub page https://github.com/sidgan/ETCI-2021-Competition-on-Flood-Detection.

Codes used to run the Paul & Ganju (2021) algorithm should be used in the following order:
1) download_s1_RTCScenes_hyp3.ipynb
2) scale_chip_s1_RTCScenes.ipynb
3) run_inference.ipynb
4) post_process_preds.ipynb

Functions used throughout the sections are stored in the folder "helpers". Note that a Google Earth Engine (GEE) account is required to run the codes for the Thomas et al. (2023) and Paul & Ganju (2021) algorithms. See: https://earthengine.google.com/.

***Many thanks and credit to Jonathan Giezendanner who created the codes to run the Thomas et al. (2023) algorithm, which were adopted from the original code by Mitchell Thomas.***

### Data
The data from this study is available via Cyverse at: https://de.cyverse.org/data/ds/iplant/home/alexsaunders/Sylhet2022Floods_local?type=folder&resourceId=1c319e32-1fff-11ee-a84e-90e2ba675364.

### Contact
For information, please contact alexsaunders@arizona.edu. 

### References

Paul, S., Ganju, S., 2021. Flood Segmentation on Sentinel-1 SAR Imagery with Semi-Supervised Learning. https://doi.org/10.48550/arXiv.2107.08369

Thomas, M., Tellman, E., Osgood, D., DeVries, B., Islam, A.S., Steckler, M.S., Goodman, M., Billah, M., 2023. A framework to assess remote sensing algorithms for satellite-based flood index insurance. IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing 1–17. https://doi.org/10.1109/JSTARS.2023.3244098
