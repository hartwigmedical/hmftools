There are pure unit tests plus more traditional functional tests in this package.

The pure unit tests use completely synthetic data.

The functional tests use a cut-down version of the Hg38 genome and a cut-down set of the v38 Ensembl data.
This data has about a dozen genes: BRAF, VHL, ZYX, etc.  They were written because it's so easy to have 
off-by-one and other errors in this domain.

Other md files in the same directory as this contain information about these genes that will
assist with understanding the functional tests.