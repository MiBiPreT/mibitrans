"""Documentation about mibitrans."""

# Add some commonly used functions as top-level imports
from mibitrans.analysis.mass_balance import mass_balance
from mibitrans.data.read import HydrologicalParameters
from mibitrans.data.read import AdsorptionParameters
from mibitrans.data.read import DegradationParameters
from mibitrans.data.read import SourceParameters
from mibitrans.data.read import ModelParameters
from mibitrans.transport.domenico import NoDecay
from mibitrans.transport.domenico import LinearDecay
from mibitrans.transport.domenico import InstantReaction
from mibitrans.visualize.plot_line import centerline
from mibitrans.visualize.plot_surface import plume_2d
from mibitrans.visualize.plot_surface import plume_3d
from mibitrans.visualize.show_mass_balance import visualize_mass_balance

import logging

logging.getLogger(__name__).addHandler(logging.NullHandler())

__author__ = "Alraune Zech"
__email__ = "a.zech@uu.nl"
__version__ = "0.4.0"
