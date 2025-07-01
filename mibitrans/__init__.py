"""Documentation about mibitrans."""
import logging

logging.getLogger(__name__).addHandler(logging.NullHandler())

__author__ = "Alraune Zech"
__email__ = "a.zech@uu.nl"
__version__ = "0.3.0"

# Add some commonly used functions as top-level imports

from mibitrans.analysis.mass_balance import MassBalance
from mibitrans.data.check_input import CheckInput
from mibitrans.data.read import from_dict
from mibitrans.transport.analytical_equation import Transport
from mibitrans.visualize.plot_line import Lineplot
from mibitrans.visualize.plot_surface import Plume
from mibitrans.visualize.show_mass_balance import generate_mass_balance_tables
from mibitrans.visualize.show_mass_balance import visualize_mass_balance




