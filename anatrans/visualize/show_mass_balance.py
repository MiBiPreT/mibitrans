"""Author: Jorrit Bakker.

Module calculating methods for mass balance visualization.
"""

import numpy as np
from prettytable import PrettyTable
from anatrans.data.parameter_information import mass_balance_renaming_dictionary as rename_dict


def generate_mass_balance_tables(mass_dict):
    table_source = PrettyTable(["Source decay", "Mass (g)"])
    table_nodecay = PrettyTable(["No decay", "Mass (g)"])
    table_lindecay = PrettyTable(["Linear decay", "Mass (g)"])
    table_instant = PrettyTable(["Instant reaction", "Mass (g)"])

    print_source = False
    print_nodecay = False
    print_lindecay = False
    print_instant = False

    for key, item in mass_dict.items():
        round_item = np.round(item, 1)
        if key != "time":#"time":
            try:
                rename_key = rename_dict[key]
            except KeyError:
                print(f"The following key is not recognized: '{key}' and is thus ignored")

            if key in ["source_mass_0", "source_mass_t", "source_mass_change"]:
                print_source = True
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
                         "plume_mass_instant_reaction", "plume_mass_degraded_instant", "electron_acceptor_mass_change"]:
                print_instant = True
                if rename_key == "source mass t = ":
                    rename_key = rename_key + str(mass_dict["time"])
                table_instant.add_row([rename_key, round_item])

    table_source = table_source if print_source else False
    table_nodecay = table_nodecay if print_nodecay else False
    table_lindecay = table_lindecay if print_lindecay else False
    table_instant = table_instant if print_instant else False

    return(table_source, table_nodecay, table_lindecay, table_instant)

def visualize_mass_balance(mass_dict):

    table_source, table_nodecay, table_lindecay, table_instant = generate_mass_balance_tables(mass_dict)

    print(f"MASS BALANCE FOR t = {mass_dict['time']}")

    if table_source is not False:
        print(table_source)
    if table_nodecay is not False:
        print(table_nodecay)
    if table_lindecay is not False:
        print(table_lindecay)
    if table_instant is not False:
        print(table_instant)

def mass_balance_excel(mass_dict):
    print("This is not yet implemented")

