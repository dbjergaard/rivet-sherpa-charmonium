. ~/rivet/local/rivetenv.sh
export PYTHIA8DATA=$(pythia8-config --xmldoc)
export RIVET_ANALYSIS_PATH=${PWD}
export LD_LIBRARY_PATH=${PWD}:${LD_LIBRARY_PATH}
