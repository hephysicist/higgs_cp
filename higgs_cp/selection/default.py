from operator import and_
from functools import reduce
from collections import defaultdict

from columnflow.production.util import attach_coffea_behavior
from columnflow.production.categories import category_ids
from columnflow.production.processes import process_ids
from columnflow.production.cms.mc_weight import mc_weight
# from columnflow.production.cms.pileup import pu_weight
# from columnflow.production.cms.scale import murmuf_weights
from columnflow.selection.cms.json_filter import json_filter
from columnflow.selection import Selector, SelectionResult, selector
from columnflow.selection.stats import increment_stats
from columnflow.columnar_util import remove_ak_column
from higgs_cp.selection.trigger import trigger_selection
from higgs_cp.selection.lepton  import study_muon_selection, study_tau_selection, mutau_selection, extra_lepton_veto, dilepton_veto
from higgs_cp.selection.jet_veto import jet_veto
from higgs_cp.production.mutau_vars import mT 

from columnflow.util import maybe_import, dev_sandbox
from higgs_cp.production.example import cutflow_features

from IPython import embed

np = maybe_import("numpy")
ak = maybe_import("awkward")


@selector(
    uses={
        "event",
        attach_coffea_behavior, 
        json_filter,
        mc_weight,
        trigger_selection,
        jet_veto,
        study_muon_selection,
        study_tau_selection,
        mutau_selection,
        extra_lepton_veto,
        dilepton_veto,
        cutflow_features,
        process_ids,
        category_ids,
        increment_stats,
        mT
    },
    produces={
        attach_coffea_behavior, 
        json_filter,
        trigger_selection,
        jet_veto,
        mc_weight,
        study_muon_selection,
        study_tau_selection,
        mutau_selection,
        extra_lepton_veto,
        dilepton_veto,
        cutflow_features,
        process_ids,
        category_ids,
        increment_stats,
        mT
    },
    sandbox=dev_sandbox("bash::$CF_BASE/sandboxes/venv_columnar_dev.sh"),
    exposed=True,
)
def default(
    self: Selector,
    events: ak.Array,
    stats: defaultdict,
    **kwargs,
) -> tuple[ak.Array, SelectionResult]:
    
    # ensure coffea behaviors are loaded
    events = self[attach_coffea_behavior](events, **kwargs)
    events = self[mT](events, **kwargs)
    
    
      # add corrected mc weights
    if self.dataset_inst.is_mc:
        events = self[mc_weight](events, **kwargs)
        
    results = SelectionResult()
    # filter bad data events according to golden lumi mask
    if self.dataset_inst.is_data:
        events, json_filter_results = self[json_filter](events, **kwargs)
        results += json_filter_results
    embed()
    events, trigger_results = self[trigger_selection](events, call_force=True, **kwargs)
    results += trigger_results
    
    events, jet_veto_results = self[jet_veto](events, call_force=True, **kwargs)
    results += jet_veto_results
    
    events, muon_results, muon_mask = self[study_muon_selection](events,
                                                                 call_force=True,
                                                                 **kwargs)
    results += muon_results
    events, tau_results, tau_mask = self[study_tau_selection](events,
                                                              call_force=True,
                                                              **kwargs)
    results += tau_results
    
    events, mutau_results, pair_mu_idx, pair_tau_idx = self[mutau_selection](events,
                                                  muon_mask,
                                                  tau_mask,
                                                  call_force=True,
                                                  **kwargs)
    results += mutau_results

    events, dilepton_veto_results = self[dilepton_veto](events,
                                                            pair_mu_idx,
                                                            pair_tau_idx,
                                                            call_force=True,
                                                            **kwargs)
    results += dilepton_veto_results
    
    events, extralep_veto_results = self[extra_lepton_veto](events,
                                                            pair_mu_idx,
                                                            pair_tau_idx,
                                                            call_force=True,
                                                            **kwargs)
    
    results += extralep_veto_results
   # write out process IDs
    events = self[process_ids](events, **kwargs)
    events = self[category_ids](events, results=results, **kwargs)
    event_sel = reduce(and_, results.steps.values())
    results.event = event_sel
    events = remove_ak_column(events, "Muon.mT")
    # some cutflow features
    events = self[cutflow_features](events, results.objects, **kwargs)
    
    # keep track of event counts, sum of weights
    # in several categories
    weight_map = {
        "num_events": Ellipsis,
        "num_events_selected": event_sel,
    }
    
    if self.dataset_inst.is_mc:
        weight_map["sum_mc_weight"] = events.mc_weight
        weight_map["sum_mc_weight_selected"] =  (events.mc_weight, event_sel)

    group_map = {}
    group_combinations = []
    
    group_map = {
        **group_map,
            # per process
            "process": {
                "values": events.process_id,
                "mask_fn": (lambda v: events.process_id == v),
            },
    }
    group_combinations.append(("process",))
            
    events, results = self[increment_stats](
        events,
        results,
        stats,
        weight_map=weight_map,
        group_map=group_map,
        group_combinations=group_combinations,
        **kwargs,
    )
    
    return events, results
