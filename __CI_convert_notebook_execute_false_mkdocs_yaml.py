import sys
import yaml

if len(sys.argv) != 2:
    sys.exit('Usage: python __CI_convert_notebook_execute_false_mkdocs_yaml.py OUTPUTFILE.yml')

outfname = sys.argv[1]
print(f"Setting mkdocs-jupyter execute flag to false in mkdocs.yml, writing output to {outfname}")


with open('mkdocs.yml', 'r') as infile:
    conf = yaml.safe_load(infile)

# If plugins block contains a dict with key 'mkdocs-jupyter', set the 'execute' flag in there to false
nb_dict_list = list(filter(lambda param: 'mkdocs-jupyter' in param, conf['plugins']))
if len(nb_dict_list) == 1:
    nb_dict_list[0]['mkdocs-jupyter']['execute'] = False

    with open(outfname, 'w') as outfile:
        yaml.safe_dump(conf, outfile, sort_keys=False)
