#!/bin/bash

set_common_vars() {

version="test"
case $1 in
    "run2" )
        config=run2_UL2018_nano_tau_v10_limited
        datasets='data_ul2018_a_single_mu,data_ul2018_b_single_mu,'`
        `'data_ul2018_c_single_mu,data_ul2018_d_single_mu'
        processes="data"
    ;;
    "run3lim")
        config="run3_2022_preEE_nano_tau_v12_limited"
        datasets='st_t_bbarq,st_tbar_bq,'`
        `'st_t_wminus_to_lnu2q,st_t_wminus_to_2l2nu,'`
        `'st_tbar_wplus_to_lnu2q,st_tbar_wplus_to_2l2nu'
        processes='st'
    ;;
    "run3")
        config="run3_2022_preEE_nano_tau_v12"
        #Datasets to use
        data='data_mu_c,data_mu_d,data_mu_e,'
        bkg_ewk='wj_incl,ww,wz,zz,wj, dy_incl,'
        bkg_top='st_t_bbarq,st_tbar_bq,'`
        `'st_t_wminus_to_lnu2q,st_t_wminus_to_2l2nu,'`
        `'st_tbar_wplus_to_lnu2q,st_tbar_wplus_to_2l2nu,'
        bkg_ttbar='tt_sl,tt_dl,tt_fh'
        datasets="$data$bkg_ewk$bkg_top$bkg_ttbar"
        processes="dy_z2mumu,dy_z2tautau,vv,tt,st,wj,data"
    ;;
    *)
    echo "Unknown run argument! Choose from: [run2, run3, run3lim]"
    exit
esac
}