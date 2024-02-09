#!/bin/bash

if [ "$1" == "run2" ]; then
    config=run2_UL2018_nano_tau_v10_limited
    dataset=data_ul2018_a_single_mu
    shift
elif [ "$1" == "run3" ]; then
    config=run3_2022_postEE_nano_tau_v12_limited
    dataset=dy_incl,signal
    shift
else
    echo "You need to choose [run2, run3] as the first argument"
    exit
fi
args=(
        --config $config
        --processes h_ggf_tautau,dy_lep
        --version test
        --datasets $dataset
        --variables muon_phi,tau_phi,mutau_mass,muon_mT
        --skip-ratio
        "$@"
    )

law run cf.PlotVariables1D "${args[@]}"