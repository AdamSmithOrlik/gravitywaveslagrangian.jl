#!/bin/bash

julia -e 'using Pkg; Pkg.activate("."); using gravitywaveslagrangian; gravitywaveslagrangian.run(:run1_dictionary)' > output.txt
