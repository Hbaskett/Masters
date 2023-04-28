# Masters
Plotting scripts and changes to ATMO and VULCAN source code

ATMO source code:

chem_neq.f90 - chemical kinetics module search get_phi to find the new kinetics module

chemistry.f90 - contains coefficients for thermal and buoyancy diffusion

grid.f90 - new set of arrays for vz

opacity.f90 - have commented out species which aren't in the C-H-N-O-S network to compute transmission spectra

VULCAN source code:

op.py - search "adection" to find functions for solving Jacobian matrices, components of transport flux have been split apart so they can be written to files to be plotted

Plotting scripts:

read_ATMO_VULCAN.py - plotting script for chemical abundances calculated in ATMO and VULCAN

read_phi.py - plotting script for transport flux components calculated in ATMO and VULCAN

read_trans.py - plotting script for transmission spectra adapted from ATMO source code

read_pt.py - plotting script for PT structure adapted from ATMO source code

read_nasa9.py - read NASA9 coefficients from VULCAN and write them into ATMO

write_atmo.py - write PT structure for ATMO format

util_commons.py - set of dictionaries for plotting scripts

Aele.py - change VULCAN elemental abundances to ATMO format
