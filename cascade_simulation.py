"""
This script simulates a cascade search on already existing XiSearch results.

This is done by carrying out a xifdr evaluation of the XiResults and removing spectra that satisfy FDR-conditions from
all successive XiResults PSMs. This process is repeated for every successive XiResult.
"""


