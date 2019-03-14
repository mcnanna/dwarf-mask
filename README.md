# dwarf-mask
Masking for dwarf searches in PS1 and DES

## Making the mask
The code for this is healpix_mask.py. The mask has \n
NSIDE = 4096
COOORD = CEL
ORDERING = NESTED

It creates a bit mask, with 0 being a good pixel and non-zero being bad for some reason. The meaning of each bit is likely to change, but right now:
0001: E(B-V) > 0.2
0010: Within 0.1 degree of object in known object catalogs
0100: Within 0.1 degree of a star in the Yale BSC

## Applying the mask
mask_check.py uses the mask to filter candidate_list.fits, and make some diagnostic plots. The purpose of this code is mostly to check the mask was made correctly and to tweak the mask, but it can also be used as a template for how it can be applied. 
