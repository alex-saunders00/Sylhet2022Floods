# Sylhet2022Floods
This repository contains code used in the following paper, which was submitted to the Institute of Electrical and Electronics Engineers (IEEE) International Symposium on Geoscience and Remote Sensing (IGARSS) 2023 conference.

## A comparison of remote sensing approaches to assess the devastating May-June 2022 floods in Sylhet, Bangladesh

Saunders, A.1, Giezendanner, J.1, Tellman, B.1, Islam, A.1, Bhuyan, A.2, Islam, A.K.M.S.3

1 Social Pixel Lab, School of Geography, Development and Environment, University of Arizona, USA

2 Flood Forecasting and Warning Centre, Bangladesh Water Development Board

3 Institute of Water and Flood Management, Bangladesh University of Engineering and Technology

## Description of the code
The code is divided into the following sections:
* Data Preparation (folder "prep")
* Data Analysis (folder "analysis")
* Results Analysis (folder "results")

The Data Prepration folder is split into a subsequent two folders for creating surface water maps using (1) the Thomas et al. (2023) ("local", non-machine learning) and (2) the Paul & Ganju (2021) ("global", machine learning) algorithms. Code for the Paul & Ganju (2021) algorithm was adapted from the GitHub page associated to the paper, which is provided in the references below.

Functions used throughout the sections are stored in the folder "helpers".

Many thanks and credit to Jonathan Giezendanner who created the codes to run the Thomas et al. (2023) algorithm, which were adopted from the original code by Mitchell Thomas.

Note that a Google Earth Engine (GEE) account is required to run the Thomas et al. (2023) algorithm. See: https://earthengine.google.com/.


