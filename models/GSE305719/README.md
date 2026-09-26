# GSE305719 reconstruction status

The project archive contains the complete inputs and code required to
generate the GSE305719 P5-P45 CORDA reconstruction family, but a complete
set of the nine frozen JSON model binaries was not located.

This is recorded explicitly rather than fabricating or mislabelling
model files.

Tracked reproducibility inputs include:

- public GSE305719 raw count matrix;
- normalized GSE305719 expression matrix;
- GSE305719 sample annotation;
- MyGene mapping cache when available;
- iMM1865 parent GEM;
- percentile/confidence reconstruction code;
- CORDA reconstruction code;
- solver/environment provenance;
- final frozen downstream result/source tables used by the thesis;
- final figure-generation code and figure source tables.

The P5-P45 models are therefore generated artifacts of the reconstruction
workflow. Full CORDA replay requires the documented Gurobi environment
and an independently supplied valid Gurobi licence. No licence file is
stored in the repository.

Published thesis endpoint sizes are P5 = 7926 reactions and
P45 = 6299 reactions and can be used as replay validation targets.
