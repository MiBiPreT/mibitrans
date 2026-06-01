#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Created on Fri Apr 10 11:33:06 2026.

@author: alraune
"""
import importlib.resources
import json
import numpy as np


class BioscreenData:
    """Class loading .json file with example output data from BIOSCREEN to compare with mibitrans models."""
    def __init__(self):
        """Initialize the example data arrays. Corresponding with model input specified in benchmarking_BIOSCREEN.ipynb.

        Properties:
            nodecay (np.ndarray) : Three-dimensional array with concentrations for no decay model,
                indexed as C(t,y,x).
            lineardecay (np.ndarray) : Three-dimensional array with concentrations for linear decay model,
                indexed as C(t,y,x).
            instant (np.ndarray) : Three-dimensional array with concentrations for instant reaction model,
                indexed as C(t,y,x).

        """
        with importlib.resources.open_text("mibitrans.data", "example_data.json") as data:
            example_data_bioscreen = json.load(data)["bioscreen"]
            self.nodecay = np.array(example_data_bioscreen["nodecay"])
            self.lineardecay = np.array(example_data_bioscreen["lineardecay"])
            self.instant = np.array(example_data_bioscreen["instant"])


class BioscreenATData:
    """Class loading .json file with example output data from BIOSCREEN-AT to compare with mibitrans models."""
    def __init__(self):
        """Initialize the example data arrays. Corresponding with model input specified in benchmarking_BIOSCREEN.ipynb.

        Properties:
            nodecay (np.ndarray) : Three-dimensional array with concentrations for no decay model,
                indexed as C(t,y,x).
            lineardecay (np.ndarray) : Three-dimensional array with concentrations for linear decay model,
                indexed as C(t,y,x).

        """
        with importlib.resources.open_text("mibitrans.data", "example_data.json") as data:
            example_data_bioscreenat = json.load(data)["bioscreenat"]
            self.nodecay = np.array(example_data_bioscreenat["nodecay"])
            self.lineardecay = np.array(example_data_bioscreenat["lineardecay"])
