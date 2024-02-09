#!/bin/bash



if [ "$1" == "run2" ]; then
    config=run2_UL2018_nano_tau_v10_limited
    dataset=data_ul2018_a_single_mu
    shift
elif [ "$1" == "run3" ]; then
    config=run3_2022_postEE_nano_tau_v12_limited
    dataset="signal"
    shift
else
    echo "You need to choose [run2, run3] as the first argument"
    exit
fi
args=(
        --config $config
        --version test
        --dataset $dataset
        "$@"
    )

law run cf.ReduceEvents "${args[@]}"