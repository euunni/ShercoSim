# DRC cosmic : 27 by 27 geometry
Repository for GEANT4 simulation &amp; analysis of the dual-readout calorimeter for general purpose.
This package is executable on CentOS7.

## How-to
### Compile
After fetching the repository, do
    
    cd install
    source envset.sh
    cd ../build
    cmake ../ -DCMAKE_INSTALL_PREFIX=../install
    make -j4 install

### Simulation

    ./bin/DRsim <macro_file> <seed> <outputfile_name>
    
### Analysis

    ./bin/analysis <path_to_root_files> <low_edge_of_hist> <truth_E> <high_edge> <outputfile_name>

e.g.)

    ./bin/analysis /home/USER/20GeV_ele_data 0 20 25 ./20GeV_ele
