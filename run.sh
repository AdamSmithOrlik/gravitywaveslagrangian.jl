#!/bin/bash

julia -e 'using Pkg; Pkg.activate("."); using gravitywaveslagrangian; gravitywaveslagrangian.run(:defualt_dictionary)' > output.txt
