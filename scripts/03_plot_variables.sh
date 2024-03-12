#!/bin/bash
source ./common.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
args=(
        --config $config
        --processes $processes
        --version $version
        --datasets $datasets
        --variables mutau_mass #muon_pt,muon_eta,muon_phi,tau_pt,tau_eta,tau_phi,muon_mT,mutau_mass,met_pt,met_phi
        --general-settings "cms-label=pw"
        "${@:2}"
    )
#echo law run cf.PlotVariables1D "${args[@]}"
law run cf.PlotVariables1D "${args[@]}"
