# sequoiaExtra

This is a collection of R functions and Fortran programs loosely associated with R package sequoia and/or for general population genetics purposes. 

Feedback and feature requests are very welcome! Or if you'd like to share some of
your own sequoia-related functions here, please do contact me!

> [!CAUTION]
> These functions and programs are in various stages of development, use at your own risk. 

# R functions
These are (currently) not included in the package because their use cases are too
rare, or because they are still too experimental. Use with caution, as they do not include input checks yet. 

To use:
- click on the file
- click on 'raw'
- copy the URL
- `source("URL")`
 
# Fortran programs
These are aimed at very large datasets, and are complementary to the sequoia stand-alone Fortran program. 

The most well-developped are:

- `pedigree_checker` : Checks if a pedigree is consistent with the genetic SNP data, by calculating for each offspring - parent - parent trio various probabilities that either or both parents truly are parents, or otherwise related.
- `find_pairs` : Finds in the genetic SNP data all sample pairs which are likely to be duplicates, or parent-offspring, or have a specified other close relationship
- `grm_tool` : Calculate various summary statistics from (potentially huge) GRMs and/or filter out pairs with R values above or below a specified threshold

# Other
Each of the Fortran programs takes genetic SNP data in the same input format as the stand-alone sequoia ( https://github.com/JiscaH/sequoia_notR ). Reformatting from PLINK ped/map or bim/bed/fam formats can be done with the 
format_SNP_data_for_sequoia.sh script. 



