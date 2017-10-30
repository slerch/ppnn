# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 11:50:02 2017

@author: sebastian
"""


## retrieve second set of auxiliary ECMWF forecast data, 
## to augment temperature forecasts and previously downloaded auxiliary data

## based on example from
## https://software.ecmwf.int/wiki/display/WEBAPI/TIGGE+retrieval+efficiency

## surface variables: (see https://software.ecmwf.int/wiki/display/TIGGE/Parameters)
# cloud cover, surface pressure, CAPE

# ECMWF forecasts from TIGGE data set: 
#   variables: 146/147/165/166/168/176/177/228039
#   all available full years, 2007-2016
#   init time 00 UTC
#   36/48 h ahead forecasts (= valid at 12 UTC and 00 UTC)
#   0.5° resolution
#   area: -10E, 30E; 30N, 70N (large part of Europe centered around Germany)
    
#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()
    
def retrieve_tigge_data():
    date1 = [str(i) + "-01-01" for i in xrange(2016,2017)] 
    date2 = [str(i) + "-12-31" for i in xrange(2016,2017)]   
    dates = date1
    for j in range(0,1):
        dates[j] = date1[j] + "/to/" + date2[j]
    data_dir = "/media/sebastian/Elements/Postproc_NN/data/forecasts/auxiliary/" 
    for date in dates:
        target = data_dir + "ecmwf_aux_surface_more_" + date[:4] + ".grib"
        tigge_request(date, target)
          
def tigge_request(date, target):
    '''
       A TIGGE request for ECMWF perturbed forecasts of auxiliary surface variables.
    '''
    server.retrieve({
        'origin'    : "ecmf",
        'levtype'   : "sfc",
        'number'    : mem_numbers,
        'expver'    : "prod",
        'dataset'   : "tigge",
        'step'      : "36/48",
        'grid'      : "0.5/0.5",
        'param'     : "146/147/165/166/168/176/177/228039",
        'area'      : "70/-10/30/30",
        'time'      : "00",
        'date'      : date,
        'type'      : "pf",
        'class'     : "ti",
        'target'    : target,
    })
 
if __name__ == '__main__':
    mem_numbers = ''.join([''.join([str(i) + "/" for i in xrange(1,50)]),'50']) 
    retrieve_tigge_data()