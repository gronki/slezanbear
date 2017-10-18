#!/bin/bash

wget --continue https://data.ngdc.noaa.gov/instruments/remote-sensing/passive/spectrometers-radiometers/imaging/viirs/dnb_composites/v10//2015/SVDNB_npp_20150101-20151231_75N060W_v10_c201701311200.tgz
wget --continue http://published-files.eea.europa.eu/eudem/entr_r_4258_1_arcsec_gsgrda-eudem-dem-europe_2012_rev1/eudem_tiles_5deg/eudem_dem_5deg_n{40,45,50,55}e0{10,15,20}.tif
