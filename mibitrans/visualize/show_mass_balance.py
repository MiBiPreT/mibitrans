"""Author: Jorrit Bakker.

Module calculating methods for mass balance visualization.
"""

import numpy as np
from prettytable import PrettyTable
from mibitrans.data.parameter_information import mass_balance_renaming_dictionary as rename_dict


def generate_mass_balance_tables(mass_dict):
    """Generate mass balance tables for relevant model components based on input dictionary.

    Args :
        mass_dict (dict) : Dictionary with mass balance with the structure as generated by the mass_balance module.

    """
    # Generate table object with object headers. Seperate tables for each component for clarity
    # table_source will contain the source decay masses from no decay and linear decay model mode.
    table_source = PrettyTable(["Source decay", "Mass (g)"])
    # table_nodecay will contain the plume masses for no decay model mode.
    table_nodecay = PrettyTable(["No decay", "Mass (g)"])
    # table_lindecay will contain the plume masses for linear decay model mode.
    table_lindecay = PrettyTable(["Linear decay", "Mass (g)"])
    # table_instant will contain the plume and source masses for instant reaction model mode.
    table_instant = PrettyTable(["Instant reaction", "Mass (g)"])
    # table_electron will contain the change in electron acceptor and byproduct masses for instant reaction model mode.
    table_electron = PrettyTable(["Electron acceptors/byproducts", "O2", "NO3-", "Fe2+", "SO4 2-", "CH4"])

    # Initialize flags which keep track of which model components are present
    print_source = False
    print_nodecay = False
    print_lindecay = False
    print_instant = False
    print_electron = False

    # Loop to add each mass balance components to their respective table
    for key, item in mass_dict.items():
        # Round mass balance to one decimal to avoid clutter
        round_item = np.round(item, 1)

        # Time does not need to be renamed and does not need to be added as a row to a table
        if key != "time":
            # Check if mass balance component is recognized, ignore if it is not.
            try:
                # Rename keys to comprehensive table headers
                rename_key = rename_dict[key]
            except KeyError:
                print(f"The following key is not recognized: '{key}' and is thus ignored")
                rename_key = None

            # Add the right mass balance model components to each mass balance table
            if key in ["source_mass_0", "source_mass_t", "source_mass_change"]:
                print_source = True
                # Add time information to name of source mass at t = t
                if rename_key == "mass t = ":
                    rename_key = rename_key + str(mass_dict["time"])
                table_source.add_row([rename_key, round_item])

            elif key in ["plume_mass_no_decay", "transport_outside_extent_nodecay"]:
                print_nodecay = True
                table_nodecay.add_row([rename_key, round_item])

            elif key in ["plume_mass_linear_decay", "transport_outside_extent_lineardecay",
                         "plume_mass_degraded_linear"]:
                print_lindecay = True
                table_lindecay.add_row([rename_key, round_item])

            elif key in ["source_mass_instant_t", "source_mass_instant_change", "plume_mass_no_decay_instant_reaction",
                         "plume_mass_instant_reaction", "plume_mass_degraded_instant"]:
                print_instant = True
                if rename_key == "source mass t = ":
                    rename_key = rename_key + str(mass_dict["time"])
                table_instant.add_row([rename_key, round_item])

            elif key == "electron_acceptor_mass_change":
                print_electron = True
                oxy, no, fe, so, ch = round_item
                # Oxygen, Nitrate and Sulfate are electron acceptors and thus consumed (negative change),
                # Iron2+ and Methane are byproducts from electron acceptors and thus generated (postive change).
                table_electron.add_row([rename_key, -oxy, -no, f"+{fe}", -so, f"+{ch}"])

    # If mass_balance components were not present in mass balance dictionary, table get set to false,
    # preventing returning tables with only headers.
    table_source = table_source if print_source else False
    table_nodecay = table_nodecay if print_nodecay else False
    table_lindecay = table_lindecay if print_lindecay else False
    table_instant = table_instant if print_instant else False
    table_electron = table_electron if print_electron else False

    return table_source, table_nodecay, table_lindecay, table_instant, table_electron

def visualize_mass_balance(mass_dict) -> None:
    """Takes dictionary with mass balance and prints it as stylized tables.

    Args :
        mass_dict (dict) : Dictionary with mass balance with the structure as generated by the mass_balance module.

    """
    # Mass balance table generating function separate to allow user to choose whether to print tables
    # or just get their objects.
    table_source, table_nodecay, table_lindecay, table_instant, table_electron = generate_mass_balance_tables(mass_dict)

    print(f"MASS BALANCE FOR t = {mass_dict['time']}")

    # If a table contained no entries, it is set to false by generating function, allowing this function to only print
    # filled tables.
    if table_source is not False:
        print(table_source)
    if table_nodecay is not False:
        print(table_nodecay)
    if table_lindecay is not False:
        print(table_lindecay)
    if table_instant is not False:
        print(table_instant)
    if table_electron is not False:
        print(table_electron)
