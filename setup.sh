. ~/rivet/local/rivetenv.sh
export PYTHIA8DATA=$(pythia8-config --xmldoc)
# hard code this so we can source it from anywhere and get the right paths
export RIVET_ANALYSIS_PATH=${HOME}/rivet/Analysis/sherpa-dev/
export LD_LIBRARY_PATH=${HOME}/rivet/Analysis/sherpa-dev/:${LD_LIBRARY_PATH}
