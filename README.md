# ArchiJulia
Translation of ArchiSimple from cpp to Julia

This is a multiprocess version of ArchiSimple ( https://www.quantitative-plant.org/model/ArchiSimple ) translated in Julia language for easier manipulation and interoperability.

The project contains a R script file to call the function with an explanation of how to use it. This version allows for looping through up to 5 parameters, without a limitation on the number of values tested.

Tested on a system with a Ryzen 7 5800X3D 8-core 16-threads processor and 64GB of ram, it takes ~35 seconds to generate 1000 independent root systems, with ~11 seconds overhead (attribution of CPUs and starting the 16 workers). Ram usage scales with the number of workers, ~3go with 8 workers and ~6go with 16.
