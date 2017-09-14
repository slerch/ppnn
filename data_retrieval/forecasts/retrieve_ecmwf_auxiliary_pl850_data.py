# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 16:27:28 2017

@author: sebastian
"""

## retrieve auxiliary ECMWF forecast data, to augment temperature forecasts
## based on example from
## https://software.ecmwf.int/wiki/display/WEBAPI/TIGGE+retrieval+efficiency

## variables at pressure level 850 hPa: 
# (see https://software.ecmwf.int/wiki/display/TIGGE/Parameters)
# u,v wind components, specific humidity

# ECMWF forecasts from TIGGE data set: 
#   all available full years, 2007-2016
#   init time 00 UTC
#   36/48 h ahead forecasts (= valid at 12 UTC and 00 UTC)
#   0.5° resolution
#   area: -10E, 30E; 30N, 70N (large part of Europe centered around Germany)
    
#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()
    
def retrieve_tigge_data():
    date1 = [str(i) + "-01-01" for i in xrange(2007,2017)] 
    date2 = [str(i) + "-12-31" for i in xrange(2007,2017)]   
    dates = date1
    for j in range(0,10):
        dates[j] = date1[j] + "/to/" + date2[j]
    data_dir = "/media/sebastian/Elements/Postproc_NN/data/forecasts/auxiliary/" 
    for date in dates:
        target = data_dir + "ecmwf_aux_pl850_" + date[:4] + ".grib"
        tigge_request(date, target)
          
def tigge_request(date, target):
    '''
       A TIGGE request for ECMWF perturbed forecasts of auxiliary variables at pressure level 850 hPa.
    '''
    server.retrieve({
        'origin'    : "ecmf",
        'levelist'  : "850",
        'levtype'   : "pl",
        'number'    : mem_numbers,
        'expver'    : "prod",
        'dataset'   : "tigge",
        'step'      : "36/48",
        'grid'      : "0.5/0.5",
        'param'     : "131/132/133",
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
