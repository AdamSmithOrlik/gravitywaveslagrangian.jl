#!/bin/bash

julia -e 'using Pkg; Pkg.activate("."); using gravitywaveslagrangian; gravitywaveslagrangian.run()' > output.txt
