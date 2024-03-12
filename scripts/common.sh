#!/bin/bash

set_common_vars() {

version="preEE"
case $1 in
    "run2" )
        config=run2_UL2018_nano_tau_v10_limited
        datasets='data_ul2018_a_single_mu,data_ul2018_b_single_mu,'`
        `'data_ul2018_c_single_mu,data_ul2018_d_single_mu'
        processes="data"
    ;;
    "run3lim")
        config="run3_2022_preEE_nano_tau_v12_limited"
        datasets='wj_incl,dy_incl,data_mu_c'
        processes='wj,dy_lep,data'
    ;;
    "run3")
        config="run3_2022_preEE_nano_tau_v12"
        datasets='data_mu_c,data_mu_d,data_mu_e,'`
        `'wj_incl,ww,wz,zz,wj,dy_incl,'`
        `'tt_sl,tt_dl,tt_fh'
        processes="data,dy_lep,ww,wz,zz,tt_sl,tt_dl,tt_fh,wj"
    ;;
    *)
    echo "Unknown run argument! Choose from: [run2, run3, run3lim]"
    exit
esac
}