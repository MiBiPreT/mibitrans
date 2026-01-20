"""Documentation about mibitrans."""

# Add some commonly used functions as top-level imports
from mibitrans.data.parameters import HydrologicalParameters
from mibitrans.data.parameters import AttenuationParameters
from mibitrans.data.parameters import SourceParameters
from mibitrans.data.parameters import ModelParameters
from mibitrans.data.parameter_information import ElectronAcceptors
from mibitrans.data.parameter_information import UtilizationFactor
from mibitrans.transport.models import Anatrans
from mibitrans.transport.models import Bioscreen
from mibitrans.transport.models import Mibitrans
from mibitrans.visualize.plot_line import centerline
from mibitrans.visualize.plot_line import transverse
from mibitrans.visualize.plot_line import breakthrough
from mibitrans.visualize.plot_surface import plume_2d
from mibitrans.visualize.plot_surface import plume_3d

import logging

logging.getLogger(__name__).addHandler(logging.NullHandler())

__author__ = "Alraune Zech"
__email__ = "a.zech@uu.nl"
__version__ = "0.5.1"
