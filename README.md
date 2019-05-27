This is hopelessly out of date. I'll update it eventually. In the meantime, just concact me if you want to know how to use the code.

# dwarf-mask
Masking and candidate filtering for dwarf searches in PS1 and DES. 
Creates two masks, one for DES and one for PS1.
Applies the masks and filters candidates coming from the output of from either the simple or ugali search algorithm. 

## Making the mask
The code for this is healpix_mask.py. The mask has
NSIDE = 4096,
COOORD = CEL,
ORDERING = NESTED.

It creates a bit mask, with 0 being a good pixel and non-zero being bad for some reason. The meaning of each bit may change, but right now:
00001: E(B-V) > 0.2
00010: Within the radius of an object in most catalogs in $UGALIDIR/catalogs
00100: Within the radius of an object in McConnachie15 and ExtraDwarf catalofs in $UGALIDIR/catalogs
01000: Within 0.1 deg of star in the Yale bright star catalog
10000: Outside of footprint

If a radius is not given in the object catalogs, it defaults to 0.1 deg. 
The footprint for Pan-STARRS is defined as dec > -25. 

## Applying the mask and filtering candidates
filter_candidates.py uses the mask plus additional cuts to filter candidate .fits files located in the "candidates" directory. It expects the files to be of the form candidate_list_survey_alg.fits, where survey is one of des or ps1, and alg is one of ugali or simple. The code makes some diagnostic plots and writes the results to .txt files, which are converted into LaTeX tables. 
In particular, it writes a signal.txt file showing how well known satellites are recovered, and a remains.txt file listing unassociated candidate hotspots. If possible, it will also write a remains_both.txt file showing hotspots which were found in both search algorithms. 
It generates a cumulative histogram of the results as a function of significance, showing counts for each level of cuts. 
