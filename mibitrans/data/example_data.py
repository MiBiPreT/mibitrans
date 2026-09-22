#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Created on Fri Apr 10 11:33:06 2026.

@author: alraune, Jorrit Bakker
"""

import importlib.resources
import json
import numpy as np

file_name = ["mibitrans.data", "example_data.json"]


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
        with importlib.resources.open_text(*file_name) as data:
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
        with importlib.resources.open_text(*file_name) as data:
            example_data_bioscreenat = json.load(data)["bioscreenat"]
            self.nodecay = np.array(example_data_bioscreenat["nodecay"])
            self.lineardecay = np.array(example_data_bioscreenat["lineardecay"])


class BiochlorData:
    """Class loading .json file with example output data from BIOCHLOR to compare with mibitrans models."""

    def __init__(self):
        """Initialize the example data arrays. Corresponding with model input specified in benchmarking_BIOCHLOR.ipynb.

        Properties:
            nodecay (np.ndarray) : Three-dimensional array with concentrations for no decay model,
                indexed as C(t,y,x).
            lineardecay (np.ndarray) : Three-dimensional array with concentrations for linear decay model,
                indexed as C(t,y,x).

        """
        with importlib.resources.open_text(*file_name) as data:
            example_data_biochlor = json.load(data)["biochlor"]
            self.pce = np.array(example_data_biochlor["pce"])
            self.tce = np.array(example_data_biochlor["tce"])
            self.dce = np.array(example_data_biochlor["dce"])
            self.vc = np.array(example_data_biochlor["vc"])
            self.eth = np.array(example_data_biochlor["eth"])
            self.t = np.array(example_data_biochlor["t"])
            self.y = np.array(example_data_biochlor["y"])
            self.x = np.array(example_data_biochlor["x"])
