# -*- coding: utf-8 -*-
"""
Created on Wed Aug 23 16:08:50 2017

@author: sebastian
"""
## requires installation of ECMWF web API software, see https://software.ecmwf.int/wiki/display/WEBAPI/ECMWF+Web+API+Home
## based on https://software.ecmwf.int/wiki/display/WEBAPI/Python+TIGGE+examples

#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
    
server = ECMWFDataServer()

server.retrieve({
    'origin'    : "ecmf",
    'levtype'   : "sfc",
    'number'    : "1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20/21/22/23/24/25/26/27/28/29/30/31/32/33/34/35/36/37/38/39/40/41/42/43/44/45/46/47/48/49/50",
    'expver'    : "prod",
    'dataset'   : "tigge",
    'step'      : "36/48",
    'grid'      : "0.5/0.5",
    'param'     : "167",
    'area'      : "70/-10/30/30",
    'time'      : "00",
    'date'      : "2016-01-01/to/2016-12-31",
    'type'      : "pf",
    'class'     : "ti",
})