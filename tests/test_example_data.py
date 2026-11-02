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

testing_massbalance_lineardecay_mbt = {
    "t": np.array([365, 730, 1095, 1460, 1825]),
    "_plume_mass_t": np.array([1628.51108295, 1963.92417588, 2016.88366043, 2018.72187961, 2016.29948848]),
    "_source_mass_t": np.array(
        [1997376.84643129, 1994757.1333299, 1992140.85618339, 1989528.01048526, 1986918.59173488]
    ),
    "_delta_source_t": np.array([2623.15356871, 5242.8666701, 7859.14381661, 10471.98951474, 13081.40826512]),
    "_degraded_mass_t": np.array([1219.80171108, 3334.71034189, 4767.88296906, 5212.29452639, 5294.24139782]),
    "_electron_acceptor_change_t": None,
    "_instant_reaction_degraded_mass_t": None,
    "source_mass_finite": True,
    "model_degradation": True,
    "model_instant_reaction": False,
}

testing_massbalance_instant_bio = {
    "t": np.array([365, 730, 1095, 1460, 1825]),
    "_plume_mass_t": np.array([871.9770543, 2313.13501054, 3510.74596612, 4295.4345529, 4482.2477845]),
    "_source_mass_t": np.array(
        [1982333.73934155, 1964823.52706592, 1947467.98477742, 1930265.74625588, 1913215.45734916]
    ),
    "_delta_source_t": np.array([17666.26065845, 35176.47293408, 52532.01522258, 69734.25374412, 86784.54265084]),
    "_degraded_mass_t": np.array([1675.89440337, 2678.48256742, 3102.50500089, 2895.02054688, 2837.061027]),
    "_electron_acceptor_change_t": {
        "oxygen": np.array([1748.09949717, 3278.04212717, 4231.80769676, 4500.89045658, 4529.0295963]),
        "nitrate": np.array([74.16179685, 139.06845388, 179.53123562, 190.94686785, 192.14064954]),
        "ferrous_iron": np.array([17586.94039579, 32979.0904915, 42574.55016134, 45281.68580557, 45564.78260523]),
        "sulfate": np.array([23731.77499191, 44501.90524154, 57449.99539844, 61102.99771354, 61485.00785284]),
        "methane": np.array([6992.39798869, 13112.16850867, 16927.23078704, 18003.56182631, 18116.11838521]),
    },
    "_instant_reaction_degraded_mass_t": np.array(
        [15392.52166487, 28864.10901806, 37262.29071018, 39631.64223556, 39879.41549936]
    ),
    "source_mass_finite": True,
    "model_instant_reaction": True,
    "model_degradation": True,
}

testing_massbalance_instant_bio_inf = {
    "t": np.array([365, 730, 1095, 1460, 1825]),
    "_plume_mass_t": np.array([908.01510775, 2491.34928403, 3916.1085068, 4992.16254439, 5461.6851394]),
    "_source_mass_t": float("inf"),
    "_delta_source_t": np.array([17744.74738066, 35489.49476131, 53234.24214197, 70978.98952263, 88723.73690328]),
    "_degraded_mass_t": np.array([1641.31413428, 2506.52396685, 2711.48622819, 2222.55141247, 1891.76958962]),
    "_electron_acceptor_change_t": {
        "oxygen": np.array([1751.28309003, 3288.39501843, 4254.88503639, 4537.5401989, 4579.7434846]),
        "nitrate": np.array([74.29685836, 139.50766745, 180.51027427, 192.50170541, 194.29214783]),
        "ferrous_iron": np.array([17618.96926937, 33083.24685205, 42806.72218428, 45650.40442525, 46074.99505719]),
        "sulfate": np.array([23774.99467674, 44642.45358349, 57763.28776675, 61600.54573046, 62173.48730609]),
        "methane": np.array([7005.13236011, 13153.58007371, 17019.54014556, 18150.16079558, 18318.9739384]),
    },
    "_instant_reaction_degraded_mass_t": np.array(
        [15420.55412074, 28955.2692199, 37465.49336015, 39954.35381666, 40325.96595349]
    ),
    "source_mass_finite": False,
    "model_instant_reaction": True,
    "model_degradation": True,
}
