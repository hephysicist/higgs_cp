
from __future__ import annotations

from columnflow.selection import Selector, SelectionResult, selector
from columnflow.util import maybe_import
from columnflow.columnar_util import set_ak_column, optional_column as optional
np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
#maybe_import("coffea.nanoevents.methods.nanoaod")


# @selector
# def trigger_selection(
#     self: Selector,
#     events: ak.Array,
#     **kwargs,
# ) -> tuple[ak.Array, SelectionResult]:

#     # start with an all-false mask
#     sel_trigger = np.array(np.zeros(len(events), dtype=np.bool_))
#     #from IPython import embed
#     #embed()
#     # pick events that passed one of the required triggers
#     for trigger in self.dataset_inst.x("require_triggers"):
#         #print(f"Requiring trigger: {trigger}")
#         single_fired = events.HLT[trigger]
#         sel_trigger = sel_trigger | single_fired

#     return events, SelectionResult(
#         steps={
#             "trigger": sel_trigger,
#         },
#     )

# @trigger_selection.init
# def trigger_selection_init(self: Selector) -> None:
#     # return immediately if dataset object has not been loaded yet
#     if not getattr(self, "dataset_inst", None):
#         return

#     # add HLT trigger bits to uses
#     self.uses |= {
#         f"HLT.{trigger}"
#         for trigger in self.dataset_inst.x.require_triggers
#     }
    

@selector(
    uses={
        # nano columns
        "TrigObj.id", "TrigObj.pt", "TrigObj.eta", "TrigObj.phi", "TrigObj.filterBits",
    },
    produces={
        # new columns
        "trigger_ids",
    },
    exposed=True,
)

def trigger_selection(
    self: Selector,
    events: ak.Array,
    **kwargs,
) -> tuple[ak.Array, SelectionResult]:
    """
    HLT trigger path selection.
    """
    any_fired = False
    trigger_data = []
    trigger_ids = []
    trig_obj_idx = []
    trig_obj_mask = []
    # index of TrigObj's to repeatedly convert masks to indices
    index = ak.local_index(events.TrigObj)

    for trigger in self.config_inst.x.triggers:
        # skip the trigger if it does not apply to the dataset
        if not trigger.applies_to_dataset(self.dataset_inst):
            continue
        # get bare decisions

        fired = events.HLT[trigger.hlt_field] == 1
        any_fired = any_fired | fired
        any_fired = ak.to_numpy(any_fired)

        # get trigger objects for fired events per leg
        leg_masks = []
        all_legs_match = True
        for leg in trigger.legs:
            # start with a True mask
            leg_mask = abs(events.TrigObj.id) >= 0
            # pdg id selection
            if leg.pdg_id is not None:
                leg_mask = leg_mask & (abs(events.TrigObj.id) == leg.pdg_id)
            # pt selection
            if leg.min_pt is not None:
                leg_mask = leg_mask & (events.TrigObj.pt >= leg.min_pt)
            # trigger bits match
            if leg.trigger_bits is not None:
                # OR across bits themselves, AND between all decision in the list
                for bits in leg.trigger_bits:
                    leg_mask = leg_mask & ((events.TrigObj.filterBits & bits) != 0)
                    
            leg_masks.append(index[leg_mask])
            trig_obj_idx.append(index[leg_mask])
            trig_obj_mask.append(leg_mask)
            # at least one object must match this leg
            all_legs_match = all_legs_match & ak.any(leg_mask, axis=1)

        # final trigger decision
        fired_and_all_legs_match = fired & all_legs_match

        # store all intermediate results for subsequent selectors
        trigger_data.append((trigger, fired_and_all_legs_match, leg_masks))
        # store the trigg er id
        ids = ak.where(fired_and_all_legs_match, np.float32(trigger.id), np.float32(np.nan))
        trigger_ids.append(ak.singletons(ak.nan_to_none(ids)))

    # store the fired trigger ids
    trigger_ids = ak.concatenate(trigger_ids, axis=1)
    events = set_ak_column(events, "trigger_ids", trigger_ids, value_type=np.int32)
    
    return events, SelectionResult(
        steps={
            "trigger": any_fired > 0,
        },
        aux={
            "trigger_data"  : trigger_data,
            "trig_obj_idx"  : trig_obj_idx,
            "trig_obj_mask" : trig_obj_mask, 
        },
    )


@trigger_selection.init
def trigger_selection_init(self: Selector) -> None:
    if getattr(self, "dataset_inst", None) is None:
        return
    # full used columns
    self.uses |= {
        optional(trigger.name)
        for trigger in self.config_inst.x.triggers
        if trigger.applies_to_dataset(self.dataset_inst)
    }
    


@selector(
    uses={
        # Muon nano columns
        f"Muon.{var}" for var in [
        "pt", "eta","phi", "dz", "dxy", "mediumId", "pfRelIso04_all",
        ] 
    } | {
        f"TrigObj.{var}" for var in [
            "id", "pt", "eta", "phi", "filterBits",
        ]
    },
    exposed=False,
)
def trigger_matching(
    self: Selector,
    events: ak.Array,
    pair_mu_idx,
    pair_tau_idx,
    trigger_results,
    **kwargs,
) -> ak.Array :
    empty_indices = ak.zeros_like(1 * events.event, dtype=np.int32)[:, np.newaxis][..., :0]
    #Prepare the mask for pair muons
    raw_muon_idx = ak.local_index(events.Muon, axis=1)
    
    (buf_idx, _) = ak.broadcast_arrays(ak.firsts(pair_mu_idx,axis=1)[:,np.newaxis],raw_muon_idx)
    pair_mu_mask = ak.fill_none(buf_idx == raw_muon_idx, False)
    pair_muons =  ak.with_name(ak.drop_none(ak.mask(events.Muon, pair_mu_mask),
                                            behavior=coffea.nanoevents.methods.vector.behavior),
                               "PtEtaPhiMLorentzVector")
    #Preselect TrigObj 
    trig_obj = ak.with_name(ak.drop_none(ak.mask(events.TrigObj, trigger_results.x.trig_obj_mask[0]),
                                         behavior=coffea.nanoevents.methods.vector.behavior),
                            "PtEtaPhiMLorentzVector")
    (broadcasted_pair_mu, _) = ak.broadcast_arrays(ak.firsts(pair_muons,axis=1)[:,np.newaxis],trig_obj)
    dr_mu2trig  = ak.fill_none(trig_obj.delta_r(broadcasted_pair_mu), 999)
    is_matched = ak.any(dr_mu2trig < 0.5, axis = 1)
    pair_mu_idx = ak.where(is_matched,
                           pair_mu_idx,
                           empty_indices)
    pair_tau_idx = ak.where(is_matched,
                            pair_tau_idx,
                            empty_indices)
    is_matched = np.array(is_matched, dtype=np.bool_) 
    return events, SelectionResult(
        steps = {
            "trigger_matching" : is_matched,
        },
        objects={
            "Muon": {   
                "Muon": pair_mu_idx
            },
            "Tau": {
                "Tau": pair_tau_idx
            }
        }
    )
    