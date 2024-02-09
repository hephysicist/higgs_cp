#!/bin/bash

if [ "$1" == "run2" ]; then
    config=run2_UL2018_nano_tau_v10_limited
    dataset=data_ul2018_a_single_mu
    shift
elif [ "$1" == "run3" ]; then
    config=run3_2022_postEE_nano_tau_v12_limited
    dataset=data_mu_f
    shift
else
    echo "You need to choose [run2, run3] as the first argument"
    exit
fi
args=(
        --config $config
        --processes data
        --version multi_pair
        --datasets $dataset
        --selector-steps "trigger,muon_pt_26,muon_eta_2p4,mediumID,muon_dxy_0p045,muon_dz_0p2,muon_iso_0p15,DeepTauVSjet,DeepTauVSe,DeepTauVSmu,tau_eta_2p3,tau_dz_0p2,tau_pt_20,single_pair,extra_lep_veto,dilep_veto"
        --shape-norm
        "$@"
    )

law run cf.PlotCutflow "${args[@]}"