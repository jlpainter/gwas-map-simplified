###############################################################################
##
## Update the mesh_hgtags file
##
## Filename: mesh_hgtags_update.py
##
## Description: Create HGTags from Mesh term/code using the UMLS hierarchy
##
##
## @author: Jeffery L Painter <jeff@jivecast.com>
## @modified: 2021-Jul-22
##
###############################################################################

import csv

###############################################################################
#                             Data Paths
###############################################################################

## Assuming this is run from the current src directory
data_path = "../data/"
UMLS_PATH   = data_path + "umls/"
OUTPUT_PATH = data_path + "output/"

# File name of the output
OUTPUT_FILE = OUTPUT_PATH + "mesh_hgtags.csv"

###############################################################################
# # Global constants
###############################################################################

# Chosen delimiter for CVS files
OUTPUT_DELIMITER = "\t"

# This should not be changed.
UMLS_DELIMITER = "|"

###############################################################################
#                              BEGIN PROGRAM
###############################################################################
if __name__ == '__main__':

    ########################################################################
    # # Setup our UMLS settings to point to local extract files
    ########################################################################
    UMLS = UMLS_PATH + "MRCONSO.RRF"
    UMLS_HIER = UMLS_PATH + "MRHIER.RRF"

    # # Dictionaries used for lookups and code matching with latest MeSH    
    mesh_code_to_aui = {}
    mesh_aui_to_code = {}
    mesh_aui_tree = {}
    mesh_path_to_aui = {}
    mesh_codes = {}
    mesh_terms = {}
    mesh_tty_included = [ "MH", "NM", "N1" ]
    
    ########################################################################
    print("Loading Mesh...")
    ########################################################################
    cnt = 0
    print(UMLS)
    reader = csv.reader(open(UMLS, 'r'), delimiter=UMLS_DELIMITER)
    for fields in reader:
        src = fields[11]
        if src == "MSH":
            aui = fields[7]
            tty = fields[12]
            code = fields[13]
            term = fields[14]
            lterm = term.lower().strip()
            
            # capture all mesh terms
            mesh_terms[lterm] = code
            if tty in mesh_tty_included:
                mesh_codes[code] = term
                mesh_code_to_aui[code] = aui
                mesh_aui_to_code[aui] = code
    
    print("Connect code to MeSH tree")
    reader = csv.reader(open(UMLS_HIER, 'r'), delimiter=UMLS_DELIMITER)
    for fields in reader:
        src = fields[4]
        if src == "MSH": 
            aui = fields[1]
            tree = fields[7]
            mesh_aui_tree[aui] = tree
            mesh_path_to_aui[tree] = aui
    
    # Now let's review all MeSH terms and try to find an appropriate HGTag for it
    print("Mesh codes: " + str(len(mesh_codes.keys())))
    print("Mesh terms: " + str(len(mesh_terms.keys())))
    print("Mesh trees: " + str(len(mesh_aui_tree.keys())))
    
    fout = open( OUTPUT_FILE, 'w' )
    header = "Mesh_ID" + OUTPUT_DELIMITER + "Mesh_Term" + OUTPUT_DELIMITER + "HG_Tag" + "\n"
    fout.write(header)
    for code in mesh_code_to_aui.keys():
        mesh_term = mesh_codes[code]
        aui = mesh_code_to_aui[code]
        
        hgtag = "Other"
        if aui in mesh_aui_tree.keys():
            my_tree = mesh_aui_tree[aui]
            path = my_tree.split(".")
            if (len(path) > 0 ):
                # Tree structure: I think these are the only two top trees we want
                #   >> Diseases [c]
                #   >> Psychiatry [f]
                # the rest are non-disease related
                # See: https://meshb.nlm.nih.gov/treeView for more info on tree view
                tree_top = path[0]
                if tree_top.startswith("C") or tree_top.startswith("F"):  
            
                    # try second level path
                    if len(path) >= 2:
                        npath = path[0] + "." + path[1]
                        if npath in mesh_path_to_aui.keys():
                            hg_aui = mesh_path_to_aui[npath]
                            hg_code = mesh_aui_to_code[hg_aui]
                            hgtag = mesh_codes[hg_code]
                            
                            output = code + OUTPUT_DELIMITER + mesh_term + OUTPUT_DELIMITER + hgtag + "\n"
                            fout.write(output)
                            
                        else:
                            print("Missing path: " + npath)
                    else:
                        pass
                else:
                    pass # non disease or psychiatric condition
                

    fout.flush()
    fout.close()
    print("Done!")
