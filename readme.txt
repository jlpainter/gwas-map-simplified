-------------------------------------------------------------------------------
Directory structure:
-------------------------------------------------------------------------------

 - gwas/src
    Contains (2) python files. The first is the GWAS-MeSH mapper and
    the second will update the hgtags file.
    
    simplified_gwas_mapper.py
        - note: ignore all prior versions, this includes multiple bug fixes

    mesh_hgtags_update.py
    
    The HG-tags update should work with stock Python 3.x and requires
    no additional libraries.  The GWAS-MeSH mapper is also Python 3.x 
    compliant, but does require that you have 'requests' and the Beautiful Soup library
    installed.
    
    https://pypi.org/project/beautifulsoup4/
    1. pip3 install bs4
    2. pip3 install requets
    3. pip3 install lxml (required for processing PubMed HTML docs)

 - gwas/data

       ./catalog     # this is where catalog files, stop words, etc are contained
       ./pubmed_html # this is where PubMed HTML files are downloaded and stored
       ./umls        # this contains a subset of the UMLS (MRCONSO, MRHIER) used
                     # for mapping terms to MeSH
       ./templates   # HTML templates for report generation
       ./output      # where output is generated

  
   
-------------------------------------------------------------------------------
Update Notes (release history):
-------------------------------------------------------------------------------
 
April 2020
-----------
 - Jeff updated code to increase coverage of this filoe to 97.6%
 - Jeff downloaded new GWAS catalog Mar-08-2020, entries increased by ~2700
 - We managed to map 89.8% of those (961 of 9480 unmapped)
 
June/July 2021
-----------
 - Latest GWAS catalog contains 18,426 entries (increase of 9800 from Mar-2020)
 - Total number of UNMAPPED studies: 3504 (19.0%)
 - still mapped 81%
 
The "pubmed_html" folder contains PubMed downloads for each study id.  Currently contains 4761 files. The Python code uses the bs4 library to parse these files and extract the relevant MeSH terms associated with each study as curated by the NLM/PubMed indexes.

Sep-30-2021
-----------
 - Noted that there were some missing maps that Alistair found in an old map file from July_2020
 
 - Reworked all logic (using lots of python scripts) to simplify the overall process, align
   maps that had multiple targets to the one with the highest frequency (another bug I had found)
   
 - Made process more clear in that the code looks at each study individually, looks at the mapped traits
   and splits those into an array, then looks at the overall disease trait field of the study.
   The code attempts to map each individual mapped trait to a Mesh term (you will see these now
   broken out in the end map file - a study CAN have multiple row entries, but the rows should
   be unique with regards to  study_id & mapped trait pairs.
   
   If a study is not mapped directly through a mapped trait, then we look at the overall study
   disease trait - if this is mapped, I then output a SINGLE map line where the individual mapped traits
   are joined back together as a single value, just as they appear in the original GWAS catalog. I
   do NOT split the study into multiple lines if it is mapped only at the disease trait level.
   
   Finally, if the study cannot be mapped, it gets a single line entry with the second column
   stating that the entry is "Unmapped" - this added column is a slight change in the original
   format of the final map file.  All other column headers in the output map file remain unchanged.
   
   
 - Added back in the supplementary data that was originally curated by Mark Hurle 
   which gives us the values for:
     	1. Suppl Concept	
 		2. Modifier	
 		3. Disease Background/Stratification, Phenotype/Biomarker, Treatment Response/Side_Effect or Interactions
 		
 - I merged all of the OTG/EFO maps into the primary disease/trait map file to reduce multiple
   file dependencies.  Once I got the mapping to where I wanted it, I regenerated the
   disease-trait mappings directly from the catalog map and stored it in the file referenced
   below.  This gets us to...
   
 - Mapped studies is now at 98% (only 390 of 19397 unmapped) in the 09-23-2021 release
   
   My main objective was to put ALL disease/trait maps into a manual map file and keep
   this in one centralized location. There are now essentially only (4) inputs to
   the simplified mapping script

   Input files all found in data/catalog directory:
   ================================================
       (1) gwas-catalog-v1.0.3-studies-rYYYY-MM-DD.tsv
           - The latest GWAS catalog file (can be updated by the Python program)

       (2) mesh_disease_trait_maps-20210930.csv
           - Complete list of GWAS disease traits mapped to Mesh terms
           - Any new maps should be added to this file (tab delimited)
           - Contains 3 columns:  Disease or Trait name, MeshID, MeshTerm

       (3) gwas_supplemntary_data.csv
           - Historical information originally curated by Mark Hurle
           - Contains 4 columns: 
           		1. GWAS Study Accession Number
           		2. Suppl Concept
           		3. Modifer
           		4. Disease Background/Stratification, Phenotype/Biomarker, Treatment Response/Side_Effect or Interactions

       (4) skip_words.txt - words that PubMed term matching should ignore (mostly countries/demographic background words)
-------------------------------------------------------------------------------
 
-------------------------------------------------------------------------------
Output files:
-------------------------------------------------------------------------------
    The hgtags script will generate a single file in the output
    directory labeled "mesh_hgtags.csv"


    The GWAS-MeSH mapper generates several files and adds the
    current day to each output file to make sure you can track 
    the history if needed (with the exception of 
    the pubmed_term_maps.csv file).
    
    The primary output will be "gwas_map_output_YYYY-mm-dd.csv"

    Secondary outputs include:
        1. gwas_unmapped_terms_with_counts_YYYY-mm-dd.html
        2. gwas_unmapped_terms_with_counts_YYYY-mm-dd.csv
        3. pubmed_term_maps.csv
        
    (1) and (2) contain the same entries but the HTML version
    is easier to quickly scan and see if there are any obvious
    new diseases or traits that SHOULD be manually mapped introduced
    by updates to the GWAS catalog that could not otherwise be
    auto-mapped.
    
    The pubmed_term_maps.csv file contains a quick summary of 
    GWAS study accession IDs together with the set of Mesh
    terms found through web scraping the PubMed website.
    
    This is populated for studies where a map could not be found
    using any of our methods in the Python code and should really
    only be referred to if you find a new disease or trait in the 
    HTML summary file that you believe should be mapped.
    
    The HTML output will include the study accession number if
    there is one and only one GWAS study that is linked
    to that disease or trait, otherwise it will include the total
    number of studies impacted by the missing map.  These
    are all displayed in descending order.
    
    If new studies are entered for which a pubmed link has NOT
    been downloaded, the program will also generate a shell
    script in the output directory with a series of "wget"
    commands that can be run and it will store those new entries
    in the pubmed_html folder. The python script will tell you
    in its output if there are new studies missing for which a 
    pubmed HTML download has not occurred yet.
    
    If you see new entries in the shell script, then you may
    want to re-run the python program again after you download
    the latest PubMed updates to make sure it was able to
    gather them all and update the pubmed_term_maps.csv file
    with the latest entries.
  
