If you want a better explanation, please contact me and I'll be happy to explain things more thouroughly. 

# Script overview

## `run_filters.py`
This is the main driving script. The significance thresholds are defined and set at the top of the script. Running the script will create four `Candidates` objects,
one for each survey/algorithm combination. It will then created fits files, tables, and diagnostic plots, saved into directories with those names. 

## `filter_candidates.py`
The `Candidates` object defined in this script does most of the work of filtering through the seeds returned by the search algorithms. It applies the healpix masks created by `healpix_mask.py`, as well as additional cuts. The relting object has a `data` attribute which is unfiltered seeds, as well as several `cut_cuttype` attributes which can be applied to the data for filtering. It then has several methods such as `write_signal` and `sighist` which will write the results to file and create plots. 

## `make_nice_tables.py`
This script reads in the results from using `filter_candidates.py` and puts them into nice LaTeX tables. 

## `healpix_mask.py`
This is the script which creates the healpix masks. 

# Details

## Healpix masks
The masks all use
NSIDE = 4096
COORD = CEL
ORDERING = NESTED

The bits correspond to
0000000: unmasked
0000001: E(B-V) > 0.2
0000010: Near cataloged object
0000100: Near known dwarf (McConnachie15 + ExtraDwarfs)
0001000: Near bright star in Yale BSC.
0010000: Outside of footprint
0100000: ugali processing failure
1000000: Spurious artifact (Used for PS1 data)

"Near" here specifically means within the half-light radius listed in the catalogs, but not smaller than 0.05 deg.
If there's no radius available, the default is 0.1 deg. 

When plotting the masks, the three "catalog" bits are combined, and the failure/artifact bits ard included in the footprint. 

## `Candidates` object

The `Candidates` object takes two inputs, the name of the survey (`des` or `ps1`) and the name of the algorithm (`ugali` or `simple`). 
It then reads in the candidate seeds, expeceting to find it at `candidates/candidate_list_survey_alg.fits`. 
Several "cut" attributes are defined based on the bitmasks:
`cut_ebv`
`cut_associate`
`cut_dwarf`
`cut_bsc`
`cut_footprint`
Two other cuts are defined:
`cut_modulus`
`cut_sig`

For convenience, the ebv, associate, bsc, footprint, and modulus cuts are combined into `cut_bulk`. 
Finally, all the cuts (bulk & dwarf & sig) are combined into `cut_final`.
I lied, there's one more cut. If requested (True by default), a `Candidates` object for the other algorithm will be created and 
compared with the one for the input algorithm. A `cut_cross` is defined by objects within 0.2 deg of each other passing the `cut_final`
for both algorithms.

The resulting `Candidates` object ends up with a `signal` array attribute. This is created by running through the objects in 
McConnachie15 and ExtraDwarfs and matching withing a 0.2 deg radius with the candidate seeds. No cuts are applied, but 
whether the object was cut and its bit on the bitmask are columns in the resulting array.

There are several methods for writing the results to file, creating LaTeX tables, and plotting. The `doitall()` method does them all. 
