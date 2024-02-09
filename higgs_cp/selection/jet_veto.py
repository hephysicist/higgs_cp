# coding: utf-8

"""
Lepton selection methods.
"""

from __future__ import annotations

from columnflow.selection import Selector, SelectionResult, selector
from columnflow.columnar_util import set_ak_column
from columnflow.util import DotDict, maybe_import
from IPython import embed




np = maybe_import("numpy")
ak = maybe_import("awkward")


@selector(
    uses={
        # Jet columns
        f"Jet.{var}" for var in [
        "pt", "eta", "btagDeepFlavB"
        ] 
    },
    exposed=False,
)
def jet_veto(
    self: Selector,
    events: ak.Array,
    **kwargs,
) -> tuple[ak.Array, ak.Array, ak.Array]:
    
    
    is2018 = self.config_inst.campaign.x.year == 2018
    is2022 = self.config_inst.campaign.x.year == 2022
    
    if   is2018: deepJetM_WP = 0.2783 # https://btv-wiki.docs.cern.ch/ScaleFactors/UL2018/
    elif is2022: deepJetM_WP = 0.3196 # https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer22EE/
    else       : deepJetM_WP = 1.0
    
    
    jet_selections = {
        "jet_veto_pt_20"           : events.Jet.pt > 20,
        "jet_veto_eta_2p4"         : abs(events.Jet.eta) < 2.4,
        "deepJet_veto_medium"      : events.Jet.btagDeepFlavB > deepJetM_WP
    }
    
    
    buffer_mask = abs(events.event) > 0

    for cut_name in jet_selections.keys():
        buffer_mask = buffer_mask & jet_selections[cut_name]
        
    b_jet_veto_evts = np.array( ~ ak.any(buffer_mask, axis=1), dtype=np.bool_) # Requiring no jets passing jet_selections
    return events, SelectionResult(
        steps = { "b_jet_veto": b_jet_veto_evts },
        
        objects = {
            "Jet" : {
                "Jet": ~ buffer_mask
            }
                
        }
        
    )
