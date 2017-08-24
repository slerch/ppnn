# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 10:44:53 2017

@author: sebastian
"""

## retrieve ECMWF forecast data, based on example from
## https://software.ecmwf.int/wiki/display/WEBAPI/TIGGE+retrieval+efficiency

# ECMWF forecasts from TIGGE data set: 
#   T2M fields 
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
    mem_numbers = ''.join([''.join([str(i) + "/" for i in xrange(1,50)]),'50'])    
    for date in dates:
        tigge_request(date)
          
def tigge_request(date):
    '''
       A TIGGE request for ECMWF perturbed forecasts of T2M.
    '''
    server.retrieve({
        'origin'    : "ecmf",
        'levtype'   : "sfc",
        'number'    : mem_numbers,
        'expver'    : "prod",
        'dataset'   : "tigge",
        'step'      : "36/48",
        'grid'      : "0.5/0.5",
        'param'     : "167",
        'area'      : "70/-10/30/30",
        'time'      : "00",
        'date'      : date,
        'type'      : "pf",
        'class'     : "ti",
    })
 
if __name__ == '__main__':
    retrieve_tigge_data()