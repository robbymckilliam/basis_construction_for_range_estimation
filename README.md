# basis_construction_for_range_estimation
Software and tex for paper

A. Akhlaq, R. G. McKilliam, and R. Subramanian. Basis construction for range estimation by phase unwrapping. IEEE Signal Process. Letters, 22(11):2152-2156, Nov 2015.

Data from the simulations already exists in this repository.  If you are using bash then

bash build.sh

will build the pdf using this data.  The code/ directory contains the source code.

bash runsim.sh

will run the simulations and regenerate data for the paper.  You need a working java jvm and also Scala 2.10.x installed. 

Makes use of the ScalaNumbers library (https://github.com/robbymckilliam/ScalaNumbers) for infinite precision matrix operations and the closest lattice point algorithms (aka sphere decoders) availabe in the JavaLatticeLibrary (https://github.com/robbymckilliam/JavaLatticeLibrary)
