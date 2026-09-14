"""Author: Jorrit Bakker.

File containing results from testing data for the transport model.
"""

import importlib.resources
import json
import numpy as np


class ExampleTestData:
    """Class loading .json file with concentration data to test functionality and output transport models."""

    def __init__(self):
        """Initialize the example data arrays. Corresponding with model input specified in conftest.py."""
        with importlib.resources.open_text("tests", "test_example_data.json") as data:
            example_test_data = json.load(data)
            self.nodecay_anatrans = np.array(example_test_data["nodecay_anatrans"])
            self.lineardecay_anatrans = np.array(example_test_data["lineardecay_anatrans"])
            self.instantreaction_anatrans = np.array(example_test_data["instantreaction_anatrans"])
            self.nodecay_mibitrans = np.array(example_test_data["nodecay_mibitrans"])
            self.lineardecay_mibitrans = np.array(example_test_data["lineardecay_mibitrans"])
            self.instantreaction_mibitrans = np.array(example_test_data["instantreaction_mibitrans"])
            self.nodecay_bioscreen = np.array(example_test_data["nodecay_bioscreen"])
            self.lineardecay_bioscreen = np.array(example_test_data["lineardecay_bioscreen"])
            self.instantreaction_bioscreen = np.array(example_test_data["instantreaction_bioscreen"])


testing_massbalance_nodecay_bio = {
    "t": np.array([365, 730, 1095, 1460, 1825]),
    "_plume_mass_t": np.array([2547.87145767, 4991.61757796, 6613.25096702, 7190.45509978, 7319.3088115]),
    "_source_mass_t": np.array(
        [1997376.84643129, 1994757.1333299, 1992140.85618339, 1989528.01048526, 1986918.59173488]
    ),
    "_delta_source_t": np.array([2623.15356871, 5242.8666701, 7859.14381661, 10471.98951474, 13081.40826512]),
    "_degraded_mass_t": None,
    "_electron_acceptor_change_t": None,
    "_instant_reaction_degraded_mass_t": None,
    "source_mass_finite": True,
    "model_degradation": False,
    "model_instant_reaction": False,
}

testing_massbalance_lineardecay_ana = {
    "t": np.array([365, 730, 1095, 1460, 1825]),
    "_plume_mass_t": np.array([1625.91225957, 1958.42344739, 2011.60695846, 2013.67822011, 2011.28905331]),
    "_source_mass_t": np.array(
        [1997376.84643129, 1994757.1333299, 1992140.85618339, 1989528.01048526, 1986918.59173488]
    ),
    "_delta_source_t": np.array([2623.15356871, 5242.8666701, 7859.14381661, 10471.98951474, 13081.40826512]),
    "_degraded_mass_t": np.array([1215.93229094, 3312.18208233, 4756.57543617, 5226.36475666, 5319.36988327]),
    "_electron_acceptor_change_t": None,
    "_instant_reaction_degraded_mass_t": None,
    "source_mass_finite": True,
    "model_degradation": True,
    "model_instant_reaction": False,
}

testing_massbalance_instant_mbt = {
    "t": np.array([365, 730, 1095, 1460, 1825]),
    "_plume_mass_t": np.array([1274.95518137, 2664.20593542, 3767.11945785, 4436.89399159, 4471.80347887]),
    "_source_mass_t": np.array(
        [1982333.73934155, 1964823.52706592, 1947467.98477742, 1930265.74625588, 1913215.45734916]
    ),
    "_delta_source_t": np.array([17666.26065845, 35176.47293408, 52532.01522258, 69734.25374412, 86784.54265084]),
    "_degraded_mass_t": np.array([1573.35761266, 2634.42858236, 3017.64717164, 2794.12241441, 2838.73740743]),
    "_electron_acceptor_change_t": {
        "oxygen": np.array([1936.17805868, 3472.4727154, 4337.57185311, 4529.31972161, 4541.44224634]),
        "nitrate": np.array([82.14088734, 147.31702429, 184.01819983, 192.15295789, 192.66724681]),
        "ferrous_iron": np.array([19479.12471152, 34935.1800458, 43638.60167375, 45567.70144168, 45689.66138744]),
        "sulfate": np.array([26285.08394807, 47141.44777265, 58885.82394531, 61488.94652371, 61653.51898065]),
        "methane": np.array([7744.7122347, 13889.89086158, 17350.28741246, 18117.27888645, 18165.76898537]),
    },
    "_instant_reaction_degraded_mass_t": np.array(
        [17048.60779573, 30576.1265814, 38193.57469641, 39881.97013662, 39988.71247299]
    ),
    "source_mass_finite": True,
    "model_degradation": True,
    "model_instant_reaction": True,
}

testing_massbalance_instant_mbt_inf = {
    "t": np.array([365, 730, 1095, 1460, 1825]),
    "_plume_mass_t": np.array([1323.58680137, 2860.86002677, 4194.058875, 5147.84088742, 5451.49338433]),
    "_source_mass_t": float("inf"),
    "_delta_source_t": np.array([17744.74738066, 35489.49476131, 53234.24214197, 70978.98952263, 88723.73690328]),
    "_degraded_mass_t": np.array([1526.73062587, 2445.16948716, 2606.17833265, 2107.91520193, 1893.34623264]),
    "_electron_acceptor_change_t": {
        "oxygen": np.array([1940.69129974, 3486.40514084, 4363.93209705, 4567.22770608, 4593.63418442]),
        "nitrate": np.array([82.33235817, 147.90809688, 185.13651321, 193.76117541, 194.88145025]),
        "ferrous_iron": np.array([19524.53065196, 35075.3486897, 43903.80170365, 45949.07873994, 46214.74391599]),
        "sulfate": np.array([26346.35461469, 47330.59100296, 59243.68422661, 62003.57613101, 62362.06407941]),
        "methane": np.array([7762.76519897, 13945.62056337, 17455.7283882, 18268.91082431, 18374.53673768]),
    },
    "_instant_reaction_degraded_mass_t": np.array(
        [17088.34818868, 30698.80561705, 38425.68426828, 40215.76090374, 40448.27758293]
    ),
    "source_mass_finite": False,
    "model_degradation": True,
    "model_instant_reaction": True,
}

testing_massbalance_instant_ana_inf = {
    "t": np.array([365, 730, 1095, 1460, 1825]),
    "_plume_mass_t": np.array([1307.95320597, 2843.96453853, 4179.04985145, 5156.20424127, 5504.06469115]),
    "_source_mass_t": float("inf"),
    "_delta_source_t": np.array([17744.74738066, 35489.49476131, 53234.24214197, 70978.98952263, 88723.73690328]),
    "_degraded_mass_t": np.array([1535.89343133, 2434.01066762, 2604.5437654, 2108.51531352, 1860.85042709]),
    "_electron_acceptor_change_t": {
        "oxygen": np.array([1916.73334639, 3443.78850861, 4332.8979254, 4553.38677448, 4582.8217764]),
        "nitrate": np.array([81.31596015, 146.10011855, 183.81991199, 193.17398437, 194.42274203]),
        "ferrous_iron": np.array([19283.49912124, 34646.59954116, 43591.57912829, 45809.83057965, 46105.96453829]),
        "sulfate": np.array([26021.10724793, 46752.03793506, 58822.37183576, 61815.67499904, 62215.27744925]),
        "methane": np.array([7666.93338555, 13775.15403444, 17331.59170161, 18213.54709793, 18331.28710558]),
    },
    "_instant_reaction_degraded_mass_t": np.array(
        [16877.39148017, 30323.55384448, 38152.4194111, 40093.88749791, 40353.07120307]
    ),
    "source_mass_finite": False,
    "model_degradation": True,
    "model_instant_reaction": True,
}
