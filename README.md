# GRACE-TWSA-downscaling-ml
This repo consists of a complete R based machine learning model created for downscaling GRACE TWSA into depth to groundwater level. The code is based on the work by Seyoum, Kwon, and Milewski, titled: Downscaling GRACE TWSA Data into High-Resolution Groundwater Level Anomaly Using Machine Learning-Based Models in a Glacial Aquifer System published on Remote Sensing Journal, published here: https://doi.org/10.3390/rs11070824. The original code was based in MATLAB, and only predicting Groundwater Storage Anomaly. However, here in a new and improved form, you could be able to predict depth to groundwater level based on GRACE statellite and other datset.

The Gradient Boosting Machines (gmb) algorithm was used. I would like to acknowledge to multiple sources of information avaiable online while working on this project. The same data published on the article was used as an example. The location of the study is Illinois, USA. The data consists of 14 feature variables: month, year, Discharge, precipitation (TRMM), Temperature(LST), Vegetation, Canopy water(Canopy), SoilMoisture, Snow Water Equivalent(SWE), Sediment Thickness(Thickness), Hydraulic Conductivity, Location of the well[Latitude(Lat), Longitude(Long)], GRACE TWSA (GRACE_TWSA), and label variable: Depth to Groundwater level (dgwl).       



**Contents**

This repo is organized as follows:

