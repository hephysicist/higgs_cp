#!/bin/bash
source ./common.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
args=(
        --config $config
        --processes dy_z2mumu,dy_z2tautau #$processes
        --version $version
        --datasets $datasets
        --workflow local
        --selector-steps "trigger,muon_pt_26,muon_eta_2p4,mediumID,muon_dxy_0p045,muon_dz_0p2,muon_iso_0p15,DeepTauVSjet,DeepTauVSe,DeepTauVSmu,tau_eta_2p3,tau_dz_0p2,tau_pt_20,single_pair,extra_lep_veto,dilep_veto"
        "${@:2}"
    )
law run cf.PlotCutflow "${args[@]}"