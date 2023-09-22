###############################################################################
##
## GWAS-MeSH mapping Utility
##
## Filename: simplified_gwas_mapper.py
##
## Description: GWAS mapper most of the prior code and logic have
##              been simplified now. Using map files to store
##              prior results and apply to latest GWAS catalog
##
## @author: Jeffery L Painter <jeff@jivecast.com>
##
##
###############################################################################

from bs4 import BeautifulSoup
import csv
import datetime
import re
import requests
import os

###############################################################################
# When are we running this code?
###############################################################################
dt = datetime.datetime.today()
today = str(dt.year) + "-" + str(dt.month) + "-" + str(dt.day)

###############################################################################
#                             Data Paths
###############################################################################

## Assuming this is run from the current src directory
data_path    = "../data/"

CAT_PATH     = data_path + "catalog/"
HTML_PATH    = data_path + "pubmed_html/"
UMLS_PATH    = data_path + "umls/"
OUTPUT_PATH  = data_path + "output/"
TEMPLATE_PATH = data_path + "templates/"



###############################################################################
#                             INPUT files
###############################################################################
#
# Latest catalog downloaded from GWAS - https://www.ebi.ac.uk/gwas/docs/file-downloads
# - note download "All studies" and not all associations (smaller file)
#
# If you want to update, the code can download the latest dictionary and will update
# this value to be the latest catalog filename.
#
latest_catalog = CAT_PATH + "gwas-catalog-v1.0.3-studies-r2022-09-30.tsv"


# Skip words contains entries of Mesh terms that
# we have identified to not be of import to mapping
# disease traits - mostly countries and other
# generic level terms
SKIP_WORDS = CAT_PATH + "skip_words.txt"


# Disease traits manually mapped to MeSH
manual_trait_file = CAT_PATH + "mesh_disease_traits_maps-20210930.csv"

# Supplementary data file
supp_data_file = CAT_PATH + "gwas_supplementary_data.csv"

# set to true to download latest GWAS catalogs
# if a new catalog file is found, the program will
# download and replaces the latest_catalog
# variable defined below with the latest version available
DOWNLOAD_NEW_GWAS_CATALOG = True
gwas_url = "https://www.ebi.ac.uk/gwas/api/search/downloads/studies_new"

# Do we want to add new catalog entries to our final output?
# This should in most cases always be true - used for debugging
ADD_NEW_CATALOG = True


###############################################################################
## Global constants
###############################################################################

# Chosen delimiter for CVS files
OUTPUT_DELIMITER = "\t"

# This should not be changed.
UMLS_DELIMITER = "|"

# set to true to download latest GWAS catalogs
# if a new catalog file is found, the program will
# download and replaces the latest_catalog
# variable defined below with the latest version available
DOWNLOAD_NEW_GWAS_CATALOG = True
gwas_url = "https://www.ebi.ac.uk/gwas/api/search/downloads/studies_new"

# Do we want to add new catalog entries to our final output?
# This should in most cases always be true - used for debugging
ADD_NEW_CATALOG = True


###############################################################################
#                             OUTPUT files
###############################################################################
DOWNLOAD_SCRIPT = OUTPUT_PATH + "download_pubmed_pages_" + today + ".sh"
PUBMED_TERM_MAPS = OUTPUT_PATH + "pubmed_term_maps.csv"

# Final output for the updated map file
OUTPUT_MAP_FILE = OUTPUT_PATH + "gwas_map_output_" + today + ".csv"

# Unmapped terms (with frequencies)
UNMAPPED_TERM_FILE = OUTPUT_PATH + "gwas_unmapped_terms_with_counts_" + today + ".csv"


###############################################################################
##                              Utility classes
###############################################################################
class UMLS_Atom():

    def __init__(self):
        self.src = ""
        self.cui = ""
        self.aui = ""
        self.term = ""
        self.code = ""
        self.tty = ""
        self.tree = ""
        return
    
    def toString(self, delimiter="|"):
        return( self.src + delimiter + self.cui + delimiter + self.aui + delimiter + self.term + delimiter + self.code + delimiter + self.tty ) 
###############################################################################


##############################################################################
## GWAS single row entry per study in catalog file
##############################################################################
class GwasEntry():
    
    def __init__(self):
        
        # Mapping details
        self.isMapped = False
        self.newStudy = False
        # Fields from original GWAS catalog
        self.Date_Added_To_Catalog = ""
        self.Pubmedid = ""
        self.First_Author = ""
        self.Date = ""
        self.Journal = ""
        self.Link = ""
        self.Study = ""

        # This is unique to the study        
        self.disease_trait = ""
        self.disease_trait_mesh_map = None
        
        # GWAS mapped traits are given as a comma delimited field
        self.mapped_traits = []
        self.uris = []
        
        # Convert those traits above into objects that have the EFO URI
        self.mapped_trait_objects = []

        # modifiers
        self.Suppl_Concept = ""
        self.Modifier = ""
        self.Disease_Background_Stratification_Or_Interactions = ""
        
        self.Initial_Sample_Size = ""
        self.Replication_Sample_Size = ""
        self.Platform_Snps_Passing_Qc = ""
        self.Association_Count = ""
        self.Study_Accession = ""
        self.Genotyping_Technology = ""
        return

    # Build the mapped trait objects
    def setMappedTraitObjects(self):
        seen = {}
        index = 0
        for my_trait in self.mapped_traits:
            my_efo = ""
            if len(self.uris) >= (index + 1):
                link = self.uris[index]
                my_efo = link.replace("http://www.ebi.ac.uk/efo/", "").strip()
            
            # Keep these unique!
            if my_trait not in seen.keys():
                traitObj = MappedTrait()
                traitObj.Mapped_Trait = my_trait
                traitObj.Mapped_Trait_EFO = my_efo
                self.mapped_trait_objects.append(traitObj)
                seen[my_trait] = 1
            index += 1
        return
    #
    # Compares only the primary fields we get from GWAS catalog
    #
    def isEqual(self, other):
        
        isEqual = True

        # Primary identifier
        if self.Study_Accession != other.Study_Accession:
            isEqual = False
        
        if self.Date_Added_To_Catalog != other.Date_Added_To_Catalog:
            isEqual = False
            
        if self.Pubmedid != other.Pubmedid:
            isEqual = False
        if self.First_Author != other.First_Author:
            isEqual = False
        if self.Date != other.Date:
            isEqual = False
        if self.Journal != other.Journal:
            isEqual = False
        if self.Link != other.Link:
            isEqual = False
        if self.Study != other.Study:
            isEqual = False
        if self.disease_trait != other.disease_trait:
            isEqual = False
        if self.Initial_Sample_Size != other.Initial_Sample_Size:
            isEqual = False
        if self.Replication_Sample_Size != other.Replication_Sample_Size:
            isEqual = False
        if self.Platform_Snps_Passing_Qc != other.Platform_Snps_Passing_Qc:
            isEqual = False

        # not used in the spreadsheet            
        #         if self.Association_Count != other.Association_Count:
        #             isEqual = False
        
        if self.Mapped_Trait_Uri != other.Mapped_Trait_Uri:
            isEqual = False
            
        if self.Genotyping_Technology != other.Genotyping_Technology:
            isEqual = False
            
        return isEqual

    #
    # Get output in the same order as spreadsheet
    #
    def toCSVOutput(self):
        output = ""
        
        ##############################################################################
        # Step 1: Check if this entry has a map and update if true
        ##############################################################################
        has_disease_map_only = False
        if self.disease_trait_mesh_map != None:
            self.isMapped = True
            has_disease_map_only = True
        for mt in self.mapped_trait_objects:
            if mt.mesh_atom != None:
                self.isMapped = True
                has_disease_map_only = False
        
        ##############################################################################
        # Step 2: If there is no map, generate a single row entry in the output file
        ##############################################################################
        if self.isMapped == False:
            
            output = output + str(self.Study_Accession) + OUTPUT_DELIMITER + "Unmapped" + OUTPUT_DELIMITER
    
            # Test for mesh map
            mesh_entry = self.disease_trait_mesh_map
            output = output + OUTPUT_DELIMITER
            output = output + OUTPUT_DELIMITER
            output = output + OUTPUT_DELIMITER
            
            output = output + str(self.Suppl_Concept) + OUTPUT_DELIMITER
            output = output + str(self.Modifier) + OUTPUT_DELIMITER
            output = output + str(self.Disease_Background_Stratification_Or_Interactions) + OUTPUT_DELIMITER
            output = output + str(self.Study) + OUTPUT_DELIMITER
            output = output + str(self.disease_trait) + OUTPUT_DELIMITER
            
            # Put these all on a single line like the original
            output = output + str(",".join(self.mapped_traits)) + OUTPUT_DELIMITER
            output = output + str(",".join(self.uris)) + OUTPUT_DELIMITER
            
            output = output + str(self.Date_Added_To_Catalog) + OUTPUT_DELIMITER
            output = output + str(self.Pubmedid) + OUTPUT_DELIMITER
            output = output + str(self.Genotyping_Technology) + OUTPUT_DELIMITER
            output = output + str(self.First_Author) + OUTPUT_DELIMITER
            output = output + self.Date + OUTPUT_DELIMITER
            output = output + self.Journal + OUTPUT_DELIMITER
            output = output + self.Link + OUTPUT_DELIMITER
            output = output + str(self.Initial_Sample_Size) + OUTPUT_DELIMITER
            output = output + str(self.Replication_Sample_Size) + OUTPUT_DELIMITER
            output = output + self.Platform_Snps_Passing_Qc + OUTPUT_DELIMITER
            
            # Is it a new study?
            if self.newStudy == True:
                output = output + "Yes"
            output = output + "\n"            
                            
        else:
            ##############################################################################
            # IF none of the GWAS mapped traits could be linked directly to a single
            # mesh term, we will evaluate if the study Disease_Trait has been mapped
            # and generate a single mapping in the output file
            ##############################################################################
            if has_disease_map_only == True:
                
                output = output + str(self.Study_Accession) + OUTPUT_DELIMITER + OUTPUT_DELIMITER
    
                # Test for mesh map
                mesh_entry = self.disease_trait_mesh_map
                output = output + str(mesh_entry.code) + OUTPUT_DELIMITER
                output = output + str(mesh_entry.term) + OUTPUT_DELIMITER
                output = output + str(mesh_entry.tree) + OUTPUT_DELIMITER
                
                output = output + str(self.Suppl_Concept) + OUTPUT_DELIMITER
                output = output + str(self.Modifier) + OUTPUT_DELIMITER
                output = output + str(self.Disease_Background_Stratification_Or_Interactions) + OUTPUT_DELIMITER
                output = output + str(self.Study) + OUTPUT_DELIMITER
                output = output + str(self.disease_trait) + OUTPUT_DELIMITER
                
                # Put these all on a single line like the original
                output = output + str(",".join(self.mapped_traits)) + OUTPUT_DELIMITER
                output = output + str(",".join(self.uris)) + OUTPUT_DELIMITER
                
                output = output + str(self.Date_Added_To_Catalog) + OUTPUT_DELIMITER
                output = output + str(self.Pubmedid) + OUTPUT_DELIMITER
                output = output + str(self.Genotyping_Technology) + OUTPUT_DELIMITER
                output = output + str(self.First_Author) + OUTPUT_DELIMITER
                output = output + self.Date + OUTPUT_DELIMITER
                output = output + self.Journal + OUTPUT_DELIMITER
                output = output + self.Link + OUTPUT_DELIMITER
                output = output + str(self.Initial_Sample_Size) + OUTPUT_DELIMITER
                output = output + str(self.Replication_Sample_Size) + OUTPUT_DELIMITER
                output = output + self.Platform_Snps_Passing_Qc + OUTPUT_DELIMITER
                
                # Is it a new study?
                if self.newStudy == True:
                    output = output + "Yes"
                output = output + "\n"            
                
            else:
                ##############################################################################
                # Mapped traits exist and we will create a single row per mapped trait to 
                # output.  This could result in some mapped traits not being expressed
                # in the final map file if they are not linked to a particular Mesh term
                # However, we will capture all of those in the unmapped term output file
                # for further evaluation.
                ##############################################################################
                for mt in self.mapped_trait_objects:
                    
                    mesh_entry = None
                    if mt.mesh_atom != None:
                        # Link to mesh
                        mesh_entry = mt.mesh_atom
                        output = output + str(self.Study_Accession) + OUTPUT_DELIMITER + OUTPUT_DELIMITER
            
                        # Test for mesh map
                        if mesh_entry != None:
                            output = output + str(mesh_entry.code) + OUTPUT_DELIMITER
                            output = output + str(mesh_entry.term) + OUTPUT_DELIMITER
                            output = output + str(mesh_entry.tree) + OUTPUT_DELIMITER
                        else:
                            output = output + OUTPUT_DELIMITER
                            output = output + OUTPUT_DELIMITER
                            output = output + OUTPUT_DELIMITER
                        
                        output = output + str(self.Suppl_Concept) + OUTPUT_DELIMITER
                        output = output + str(self.Modifier) + OUTPUT_DELIMITER
                        output = output + str(self.Disease_Background_Stratification_Or_Interactions) + OUTPUT_DELIMITER
                        output = output + str(self.Study) + OUTPUT_DELIMITER
                        output = output + str(self.disease_trait) + OUTPUT_DELIMITER
                        
                        # Mapped traits are objects with the trait and EFO term
                        output = output + str(mt.Mapped_Trait) + OUTPUT_DELIMITER
                        output = output + str(mt.Mapped_Trait_EFO) + OUTPUT_DELIMITER
                        
                        output = output + str(self.Date_Added_To_Catalog) + OUTPUT_DELIMITER
                        output = output + str(self.Pubmedid) + OUTPUT_DELIMITER
                        output = output + str(self.Genotyping_Technology) + OUTPUT_DELIMITER
                        output = output + str(self.First_Author) + OUTPUT_DELIMITER
                        output = output + self.Date + OUTPUT_DELIMITER
                        output = output + self.Journal + OUTPUT_DELIMITER
                        output = output + self.Link + OUTPUT_DELIMITER
                        output = output + str(self.Initial_Sample_Size) + OUTPUT_DELIMITER
                        output = output + str(self.Replication_Sample_Size) + OUTPUT_DELIMITER
                        output = output + self.Platform_Snps_Passing_Qc + OUTPUT_DELIMITER
                        
                        # Is it a new study?
                        if self.newStudy == True:
                            output = output + "Yes"
                        output = output + "\n"
                               
        return output
##############################################################################

##############################################################################
# Each GWAS mapped trait will be split into these objects so we can 
# quickly see if we have linked it to a Mesh term (UMLS atom)
##############################################################################
class MappedTrait():

    def __init__(self):
        self.mesh_atom = None
        self.Mapped_Trait = ""
        self.Mapped_Trait_EFO = ""
        return

##############################################################################
# This method loads in our manual disease trait to Mesh mappings
##############################################################################
def getMeshMap(map_file):
    terms = {}
    lc = 0
    has_header = True
    reader = csv.reader(open(map_file, 'r', encoding = 'utf8'), delimiter="\t")
    for fields in reader:
        if has_header == True:
            has_header = False
        else:
            try:
                trait = fields[0].strip().lower()
                mesh_id = fields[1]
                terms[trait] = mesh_id
            except:
                print("Error on line: " + str(lc))
        lc = lc + 1
    return terms

##############################################################################
# This method loads in the supplementary data from an external file and
# updates the GWAS studies loaded into the dictionary
##############################################################################
def loadSupplementaryData(supplement_file, gwas_studies):
    reader = csv.reader(open(supplement_file, 'r', encoding = 'utf8'), delimiter="\t")
    for fields in reader:
        study_id = fields[0].strip()
        if study_id in gwas_studies.keys():
            entry = gwas_studies[study_id]
            entry.Suppl_Concept = fields[1]
            entry.Modifier = fields[2]
            entry.Disease_Background_Stratification_Or_Interactions = fields[3]
    return


##############################################################################
## Load GWAS studies from a standard catalog file
##############################################################################
def getGwasStudies(study_file, DELIMITER):
    
    # Expected header format (this could change if GWAS catalog format changes
    # DATE ADDED TO CATALOG    PUBMEDID    FIRST AUTHOR    DATE    JOURNAL    LINK    STUDY    DISEASE/TRAIT    INITIAL SAMPLE SIZE    REPLICATION SAMPLE SIZE    PLATFORM [SNPS PASSING QC]    ASSOCIATION COUNT    MAPPED_TRAIT    MAPPED_TRAIT_URI    STUDY ACCESSION    GENOTYPING TECHNOLOGY
    studies = {}
    has_header = True
    reader = csv.reader(open(study_file, 'r', encoding = 'utf8'), delimiter=DELIMITER)
    for fields in reader:
        if has_header == True:
            has_header = False
        else:
            
            # Generate a new entry
            entry = GwasEntry()
            entry.Study_Accession = fields[14].strip()
            if entry.Study_Accession in studies.keys():
                entry = studies[entry.Study_Accession]
            else:
                entry.Genotyping_Technology = fields[15].strip()
                entry.Date_Added_To_Catalog = fields[0].strip()
                entry.Pubmedid = fields[1].strip()
                entry.First_Author = fields[2].strip()
                entry.Date = fields[3].strip()
                entry.Journal = fields[4].strip()
                entry.Link = fields[5].strip()
                entry.Study = fields[6].strip()
                entry.disease_trait = fields[7].strip()
                entry.Initial_Sample_Size = fields[8].strip()
                entry.Replication_Sample_Size = fields[9].strip()
                entry.Platform_Snps_Passing_Qc = fields[10].strip()
                
                # Split the mapped traits and convert all to lower case
                # for easier comparisons
                my_traits = fields[12].lower().strip()
                if "," in my_traits:
                    tmp = my_traits.split(",")
                    for t in tmp:
                        t = t.lower().strip()
                        if len(t) > 0:
                            entry.mapped_traits.append(t)
                else:
                    # single trait
                    entry.mapped_traits.append(my_traits)

                # Split linked URIs for mapped traits
                # these SHOULD align with the traits above
                # and have not found a case where they didn't
                # but I did not put a lot of checks into this
                my_uris = fields[13].strip()
                if "," in my_uris:
                    links = my_uris.split(",")
                    for link in links:
                        entry.uris.append(link)
                    
                # Based on the above, generate as MappedTrait objects
                entry.setMappedTraitObjects()

            # store the entry based on it's study id as the key
            studies[entry.Study_Accession] = entry
         
    return studies


###############################################################################
## Get a single UMLS atom based on a Mesh ID
###############################################################################
def getPreferredMeshTerm(mesh_id, mesh_code_to_aui, mesh_atoms):
    
    mh_atom = None
    nm_atom = None
    n1_atom = None
    
    # Update assignments
    if mesh_id in mesh_code_to_aui.keys():
        for aui in mesh_code_to_aui[mesh_id].keys():
            atom = mesh_atoms[aui]
            if atom.tty == "MH":
                mh_atom = atom
            if atom.tty == "NM":
                nm_atom = atom
            if atom.tty == "N1":
                n1_atom = atom

    # Pick atom to return
    if mh_atom != None:
        return mh_atom
    
    if nm_atom != None:
        return nm_atom

    if n1_atom != None:
        return n1_atom
    
    return None
###############################################################################
        
    

###############################################################################
## Get a single UMLS atom based on a Mesh term
###############################################################################
def getPreferredMeshTermByTerm(mesh_term, mesh_term_to_aui, mesh_atoms):
    
    # Test that we have a string
    if mesh_term != None:
        mh_atom = None
        nm_atom = None
        n1_atom = None

        # Standardize the term as we have done with the UMLS terms
        mesh_term = mesh_term.lower().strip()
        if len(mesh_term) > 0:
        
            # Update assignments
            if mesh_term in mesh_term_to_aui.keys():
                for aui in mesh_term_to_aui[mesh_term].keys():
                    atom = mesh_atoms[aui]
                    # set based on term type from UMLS
                    if atom.tty == "MH":
                        mh_atom = atom
                    if atom.tty == "NM":
                        nm_atom = atom
                    if atom.tty == "N1":
                        n1_atom = atom
        
            # Pick atom to return in the following order
            if mh_atom != None:
                return mh_atom
            
            if nm_atom != None:
                return nm_atom
        
            if n1_atom != None:
                return n1_atom
        
    return None
###############################################################################
        
###############################################################################
def getFilename_fromCd(cd):
    """
    Gets the filename from content-disposition
    """
    if not cd:
        return None
    fname = re.findall('filename=(.+)', cd)
    if len(fname) == 0:
        return None
    return fname[0]
###############################################################################


###############################################################################
#                              BEGIN PROGRAM
###############################################################################
if __name__ == '__main__':

    ###########################################################################
    # Do we want to download the latest GWAS catalog?
    ###########################################################################
    if DOWNLOAD_NEW_GWAS_CATALOG == True:
    
        print("Downloading the latest GWAS catalog file...")
        # Download and store the latest catalog
        r = requests.get(gwas_url, allow_redirects=True)
        filename = getFilename_fromCd(r.headers.get('content-disposition'))
        outfile = CAT_PATH + filename
        if os.path.exists(outfile) == False:
            print("Saving: " + filename)
            
            open(outfile, 'wb').write(r.content)
            # Update the latest_catalog value with the current file name
            latest_catalog = outfile
        else:
            print("Latest catalog already downloded: " + filename)
            latest_catalog = outfile

    ###########################################################################
    # Skip words include MeSH terms that we do not want to
    # map to use from PubMed web-scraped links
    ###########################################################################
    skipwords = {}
    fh = open(SKIP_WORDS, 'r', encoding = 'utf8')
    for line in fh.readlines():
        term = line.strip().lower()
        if len(term) > 1:
            skipwords[term] = True
    fh.close()
    print("Skip words loaded: " + str(len(skipwords.keys())))
    
    ########################################################################
    # Load the latest catalog from GWAS and find all unique traits
    ########################################################################
    gwas_studies = getGwasStudies(latest_catalog, "\t")
    total_studies = len(gwas_studies.keys())
    print("Total studies in latest GWAS catalog: " + str(total_studies))
    unique_mapped_traits = {}
    unique_disease_trait = {}
    for study_id in gwas_studies.keys():
        entry = gwas_studies[study_id]
        for trait in entry.mapped_traits:
            unique_mapped_traits[trait] = 1
        unique_disease_trait[entry.disease_trait.lower().strip()] = 1
    
    ########################################################################        
    # Summary stats about disease traits and mapped traits from the latest catalog
    ########################################################################
    print("Total disease traits: "     + str(len(unique_disease_trait.keys())))
    print("Total GWAS mapped traits: " + str(len(unique_mapped_traits.keys())))
    
    ########################################################################
    # Add the supplementary data to studies
    ########################################################################
    loadSupplementaryData(supp_data_file, gwas_studies)
    
    ########################################################################
    ## Setup our UMLS settings to point to local extract files
    ########################################################################
    UMLS = UMLS_PATH + "MRCONSO.RRF"
    UMLS_HIER = UMLS_PATH + "MRHIER.RRF"

    ## Dictionaries used for lookups and code matching with latest MeSH
    mesh_atoms = {}
    mesh_code_to_aui = {}
    mesh_term_to_aui = {}
    mesh_tty_included = [ "MH", "NM", "N1" ]
    
    ########################################################################
    print("Loading Mesh...")
    ########################################################################
    cnt = 0
    print(UMLS)
    reader = csv.reader(open(UMLS, 'r', encoding = 'utf8'), delimiter=UMLS_DELIMITER)
    for fields in reader:
        src = fields[11]
        if src == "MSH":
            tty = fields[12]
            if tty in mesh_tty_included:
                aui = fields[7]
                code = fields[13]
                term = fields[14]
                lterm = term.strip().lower()
                
                mesh_entry = UMLS_Atom()
                mesh_entry.aui = aui
                mesh_entry.tty = tty
                mesh_entry.code = code
                mesh_entry.term = term
                mesh_atoms[aui] = mesh_entry
                
                if code not in mesh_code_to_aui.keys():
                    mesh_code_to_aui[code] = {}
                mesh_code_to_aui[code][aui] = 1
                    
                if lterm not in mesh_term_to_aui.keys():
                    mesh_term_to_aui[lterm] = {}
                mesh_term_to_aui[lterm][aui] = 1
    
    print("Connecting Mesh codes to the MeSH hierarchy")
    reader = csv.reader(open(UMLS_HIER, 'r', encoding = 'utf8'), delimiter=UMLS_DELIMITER)
    for fields in reader:
        src = fields[4]
        if src == "MSH":    
            aui = fields[1]
            try:
                mesh_atoms[aui].tree = fields[7]
            except:
                pass
    
    print("Mesh codes: " + str(len(mesh_code_to_aui.keys())))
    print("Mesh terms: " + str(len(mesh_term_to_aui.keys())))
    print("Mesh atoms: " + str(len(mesh_atoms.keys())))


    ########################################################################
    # Step 1: Map GWAS disease traits to mesh
    ########################################################################
    gwas_mapped_disease_traits_to_mesh_atoms = {}
    for trait in unique_disease_trait.keys():
        atom = None
        if trait in mesh_term_to_aui.keys():
            atom = getPreferredMeshTermByTerm(trait, mesh_term_to_aui, mesh_atoms)
            if atom != None:
                gwas_mapped_disease_traits_to_mesh_atoms[trait] = atom

    print("GWAS high level disease-traits mapped directly to MeSH: " + str(len(gwas_mapped_disease_traits_to_mesh_atoms.keys())))
        
    ########################################################################
    # Step 2: Map GWAS mapped traits to mesh
    ########################################################################
    gwas_mapped_traits_to_mesh_atoms = {}
    for trait in unique_mapped_traits.keys():
        atom = None
        if trait in mesh_term_to_aui.keys():
            atom = getPreferredMeshTermByTerm(trait, mesh_term_to_aui, mesh_atoms)
            if atom != None:
                gwas_mapped_traits_to_mesh_atoms[trait] = atom
    print("GWAS traits mapped directly to MeSH: " + str(len(gwas_mapped_traits_to_mesh_atoms.keys())))    
    
    ########################################################################
    # Step 3: Look at our manual maps and update
    ########################################################################
    manually_mapped_terms = getMeshMap(manual_trait_file)
    for trait in manually_mapped_terms.keys():
        mesh_id = manually_mapped_terms[trait]
        atom = getPreferredMeshTerm(mesh_id, mesh_code_to_aui, mesh_atoms)
        if atom != None:
            # Test for GWAS high level disease
            if trait in unique_disease_trait.keys():
                if trait not in gwas_mapped_disease_traits_to_mesh_atoms.keys():
                    # add new map!
                    gwas_mapped_disease_traits_to_mesh_atoms[trait] = atom
            
            # Test for individual mapped trait
            if trait in unique_mapped_traits.keys():
                if trait not in gwas_mapped_traits_to_mesh_atoms.keys():
                    gwas_mapped_traits_to_mesh_atoms[trait] = atom
    
    ########################################################################
    # Step 4: Use these mappings and apply to the GWAS catalog entries
    ########################################################################
    # Currently, all studies are set to unmapped status
    for study_id in gwas_studies.keys():
        entry = gwas_studies[study_id]
        
        # First look at individual traits and attempt to map at that level
        for mt in entry.mapped_trait_objects:
            term = mt.Mapped_Trait.lower().strip()
            if term in gwas_mapped_traits_to_mesh_atoms.keys():
                mt.mesh_atom = gwas_mapped_traits_to_mesh_atoms[term]
                entry.isMapped = True
            else:
                # maybe in here?
                if term in gwas_mapped_disease_traits_to_mesh_atoms.keys():
                    mt.mesh_atom = gwas_mapped_disease_traits_to_mesh_atoms[term]
                    entry.isMapped = True
            
        # Second, look at the GWAS high level disease-trait
        study_disease = entry.disease_trait.lower().strip()
        if study_disease in gwas_mapped_traits_to_mesh_atoms.keys():
            entry.disease_trait_mesh_map = gwas_mapped_traits_to_mesh_atoms[study_disease]
            entry.isMapped = True
        else:
            if study_disease in gwas_mapped_disease_traits_to_mesh_atoms.keys():
                entry.disease_trait_mesh_map = gwas_mapped_disease_traits_to_mesh_atoms[study_disease]
                entry.isMapped = True


    ########################################################################
    ## For unmapped studies, we want to build a script to make sure we can download
    ## their abstracts from PubMed
    ########################################################################
    need_to_download_new_files = False
    fout = open(DOWNLOAD_SCRIPT, 'w', encoding = 'utf8')
    unmapped = {}
    for study_id in gwas_studies.keys():
        entry = gwas_studies[study_id]
        if entry.isMapped == False:
            
            # See if we have already downloaded the file before doing it again
            html_file = HTML_PATH + entry.Study_Accession + ".html"
            if os.path.exists(html_file) == True:
                # print("File already downloaded: " + html_file)
                pass
            else:
                # build download command
                html_file = html_file.replace("../data/", "../")
                cmd = "wget -O " + html_file + " --no-clobber " + str(entry.Link)
                fout.write(cmd + "\n")
                need_to_download_new_files = True
            
            unmapped[study_id] = 1
            
    # clouse download script file
    fout.flush()
    fout.close()
    print("Unmapped studies: " + str(len(unmapped.keys())))
    if need_to_download_new_files == True:
        print("*** You need to run the shell script: " + DOWNLOAD_SCRIPT + " to download updated abstracts from PubMed.")
    
    ########################################################################
    # Extract the MeSH terms manually curated to our unmapped studies
    # and create a spreadsheet for analysis to update our manual map file
    ########################################################################
    study_terms = {}
    titles = {}
    for study_id in unmapped.keys():
             
        # get the html file       
        html_file = HTML_PATH + study_id + ".html"
        
        if os.path.exists(html_file) == True:
            content = ""
            fh = open(html_file, 'r', encoding = 'utf8')
            for line in fh.readlines():
                content += line
            fh.close()
            
            # create a dictionary for this study to map terms to
            study_terms[study_id] = {}
            soup = BeautifulSoup(content, "lxml")
            
            # extract Pub-med verbatim title for this study
            try:
                title = soup.find('title').contents[0]
                title = title.replace("- PubMed - NCBI", "").strip()
                titles[study_id] = title
                # print("Title: " + title)
                                
                # Find all links
                links = soup.findAll('a')
                for link in links:
                    if "alsec" in link.attrs:
                        if "mesh" == link.attrs["alsec"]:
                            mesh_term = link.contents[0].strip()
                            mesh_term = mesh_term.replace("*", "")
                            
                            # many terms have a split, let's get rid of it
                            if "/" in mesh_term:
                                m = mesh_term.split("/")
                                mesh_term = m[0].strip()
                            
                            mlower = mesh_term.strip().lower()
                            if mlower not in skipwords.keys():
                                study_terms[study_id][mlower] = 1      
            except:
                print("Removing broken HTML file: " + html_file)
                os.remove(html_file)
        
        else:
            print("Missing file for Study: " + study_id)

    ########################################################################
    # Generate local file for manual review
    ########################################################################
    # store frequency with which the gwas_mapped_studies trait occurs
    unmapped_trait_terms = {}
    unmapped_trait_term_to_study = {}
    for study_id in gwas_studies.keys():
        entry = gwas_studies[study_id]
        if entry.isMapped == False:
            for mt in entry.mapped_trait_objects:
                term = mt.Mapped_Trait
                if term not in  unmapped_trait_terms.keys():
                    unmapped_trait_terms[term] = 0
                    # record the study accession number
                    unmapped_trait_term_to_study[term] = entry.Study_Accession
                unmapped_trait_terms[term] = unmapped_trait_terms[term] + 1
            
    print("Total disease/traits that are unmapped: " + str(len(unmapped_trait_terms.keys())) )
    
    # Sort by frequency in descending order (biggest to smallest)
    unmapped_trait_terms_sorted = (sorted(unmapped_trait_terms.items(), key=lambda item: item[1], reverse = True))

    html_data_rows = ""    
    fout = open(UNMAPPED_TERM_FILE, 'w', encoding = 'utf8')
    header = "Row" + OUTPUT_DELIMITER + "Disease_Trait_Term" + OUTPUT_DELIMITER + "Frequency" + OUTPUT_DELIMITER + "Study Accession" + OUTPUT_DELIMITER + "Measurement"
    fout.write(header + "\n")
    row_count = 1
    for entry in unmapped_trait_terms_sorted:
        unmapped_term = entry[0]
        count = entry[1]
        is_measure = ""
        if ( "measure" in unmapped_term.lower() or " volume" in unmapped_term.lower() ):
            is_measure = "True"

        study_id = ""
        if count == 1:
            study_id = unmapped_trait_term_to_study[unmapped_term]
            
        new_row = ""
        if (len(is_measure) > 0 ):
            # comment out measures for now            
            new_row = "<tr> <td> " +str(row_count) + " </td><td> " + unmapped_term + " </td><td>" + str(count) +"</td><td> "+ study_id + "</td><td> " + is_measure +" </td> </tr>"            
        else:
            new_row = "<tr class=\"w3-green\"> <td> " +str(row_count) + " </td><td> " + unmapped_term + " </td><td>" + str(count) +"</td><td> "+ study_id + "</td><td> False </td></tr>"
            
        html_data_rows = html_data_rows + new_row + "\n"

        # Generate CVS output        
        output = str(row_count) + OUTPUT_DELIMITER + unmapped_term + OUTPUT_DELIMITER + str(count) + OUTPUT_DELIMITER + study_id + OUTPUT_DELIMITER + is_measure + "\n"
        fout.write(output)
        row_count = row_count + 1
        
    # Close the output file
    fout.flush()
    fout.close()

    ########################################################################
    # Save as HTML output for easy viewing
    ########################################################################
    html_template_file = TEMPLATE_PATH + "gwas_mesh_missing_diseases_template.html"
    html_data = ""
    with open (html_template_file, "r") as myfile:
        html_data = myfile.readlines()
    html_data = "\n".join(html_data)
    
    html_data = html_data.replace("$update_time$", today)
    html_data = html_data.replace("$csv_file_name$", "./gwas_unmapped_terms_with_counts_" + today + ".csv")
    html_data = html_data.replace("$unmapped_term_rows$", html_data_rows)
    HTML_OUTPUT_FILE = OUTPUT_PATH + "gwas_unmapped_terms_with_counts_" + today + ".html"
    fout = open(HTML_OUTPUT_FILE, 'w', encoding = 'utf8')
    fout.write(html_data)
    fout.flush()
    fout.close() 
    ########################################################################    
    
    ########################################################################
    ## Store the link between pubmed mesh terms and studies that we
    ## downloaded from the online links and processed using bs4 lib
    ########################################################################
    print("Using PubMed to help with unmapped terms...")
    fout = open(PUBMED_TERM_MAPS, 'w', encoding = 'utf8')
    header = ""
    header = header + "Unique ID" + OUTPUT_DELIMITER
    header = header + "Study Title" + OUTPUT_DELIMITER
    header = header + "Study Accession" + OUTPUT_DELIMITER
    header = header + "Disease_Trait" + OUTPUT_DELIMITER
    header = header + "Mapped_Trait" + OUTPUT_DELIMITER
    header = header + "Mapped_Trait_Freq" + OUTPUT_DELIMITER
    header = header + "PubMed_Mapped_Mesh_ID" + OUTPUT_DELIMITER
    header = header + "PubMed_Mapped_Mesh_Term"
    
    fout.write(header + "\n")
    unmapped_and_no_pubmed = {}
    for study_id in gwas_studies.keys():
        entry = gwas_studies[study_id]
        if entry.isMapped == False:
            if entry.Study_Accession in study_terms.keys():
                
                mapped_terms = study_terms[entry.Study_Accession].keys()
                if len(mapped_terms) > 0:
                    fout.write("\n" + "Study " + entry.Study_Accession + "\n")
                    for term in mapped_terms:
                        try:
                            atom = getPreferredMeshTermByTerm(term, mesh_term_to_aui, mesh_atoms)
                            mesh_id = ""
                            if atom != None:
                                mesh_id = atom.code
                                
                            output = ""
                            output = str(entry.unique_id) + OUTPUT_DELIMITER
                            output = output + titles[entry.Study_Accession] + OUTPUT_DELIMITER
                            output = output + entry.Study_Accession + OUTPUT_DELIMITER
                            output = output + entry.Disease_Trait + OUTPUT_DELIMITER
                            output = output + entry.Mapped_Trait + OUTPUT_DELIMITER
                            count = 0
                            if entry.Mapped_Trait in unmapped_trait_terms.keys():
                                count = unmapped_trait_terms[entry.Mapped_Trait]
                            output = output + str(count) + OUTPUT_DELIMITER
                            output = output + str(mesh_id) + OUTPUT_DELIMITER
                            output = output + term
                            if len(output) > 0: 
                                fout.write(output + "\n")
                        except:
                            print(" >> GWAS had mapped term: " + term + " and no Mesh map exists")
                else:
                    unmapped_and_no_pubmed[entry.Study_Accession] = 1
                    
    fout.flush()
    fout.close()
    print("Studies that are not mapped and have no PubMed recommendations: " + str(len(unmapped_and_no_pubmed.keys())))
    
    ########################################################################
    ## Update new map file
    ########################################################################
    print("Save new map file!")
    unmapped = 0
    header = "STUDY ACCESSION" + OUTPUT_DELIMITER + "HAS_MAP" + OUTPUT_DELIMITER + "MeSH ID" + OUTPUT_DELIMITER + "MeSH Name" + OUTPUT_DELIMITER + "Chosen Tree Position" + OUTPUT_DELIMITER + "Suppl Concept" + OUTPUT_DELIMITER + "Modifier" + OUTPUT_DELIMITER + "Disease Background/Stratification or Interactions" + OUTPUT_DELIMITER + "STUDY" + OUTPUT_DELIMITER + "DISEASE/TRAIT" + OUTPUT_DELIMITER + "MAPPED_TRAIT" + OUTPUT_DELIMITER + "MAPPED_TRAIT_URI" + OUTPUT_DELIMITER + "DATE ADDED TO CATALOG" + OUTPUT_DELIMITER + "PUBMEDID" + OUTPUT_DELIMITER + "GENOTYPING TECHNOLOGY" + OUTPUT_DELIMITER + "FIRST AUTHOR" + OUTPUT_DELIMITER + "DATE" + OUTPUT_DELIMITER + "JOURNAL" + OUTPUT_DELIMITER + "LINK" + OUTPUT_DELIMITER + "INITIAL SAMPLE SIZE" + OUTPUT_DELIMITER + "REPLICATION SAMPLE SIZE" + OUTPUT_DELIMITER + "PLATFORM [SNPS PASSING QC]"
    fout = open(OUTPUT_MAP_FILE, 'w', encoding = 'utf8')
    fout.write(header + "\n")
    
    for study_key in gwas_studies.keys():
        study = gwas_studies[study_key]
        fout.write(study.toCSVOutput())
        if study.isMapped == False:
            unmapped += 1

    fout.flush()
    fout.close()
    print("Total studies remaining unmapped: " + str(unmapped) + " == " + str(round( float(unmapped) / float(total_studies) , 4 ) * 100.0 ) + "%")
    print("Done!")
