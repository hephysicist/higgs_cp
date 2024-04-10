"""
Wrappers for some default sets of producers.
"""

from columnflow.production import Producer, producer
from columnflow.util import maybe_import

from columnflow.production.normalization import normalization_weights
from columnflow.production.categories import category_ids
from higgs_cp.production.example import features, cutflow_features 
from higgs_cp.production.mutau_vars import dilepton_mass, mT, rel_charge
from higgs_cp.production.weights import pu_weight, muon_weight, tau_weight
from higgs_cp.production.sample_split import split_dy
from higgs_cp.calibration.tau import tau_energy_scale
ak = maybe_import("awkward")

pdgId_map = {
    1: "down",
    2: "up",
    3: "strange",
    4: "charm",
    5: "bottom",
    6: "top",
    11: "electron",
    12: "e_neutrino",
    13: "muon",
    14: "mu_neutrino",
    15: "tau",
    16: "tau_neutrino",
    21: "gluon",
    22: "photon",
    23: "Z",
    24: "W",
    25: "Higgs",
}

@producer(
    uses={
        rel_charge, category_ids, features, normalization_weights, cutflow_features, dilepton_mass, mT, pu_weight, muon_weight, tau_weight, split_dy
    },
    produces={
        rel_charge, category_ids, features, normalization_weights, cutflow_features, dilepton_mass, mT, pu_weight, muon_weight, tau_weight, split_dy
    },
)
def default(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    events = self[rel_charge](events, **kwargs)
    events = self[category_ids](events, **kwargs)
    if self.dataset_inst.is_mc:
        events = self[normalization_weights](events, **kwargs)
        # add corrected mc weights
        processes = self.dataset_inst.processes.names()
        if ak.any(['dy' in proc for proc in processes]):
            print("Splitting Drell-Yan dataset...")
            events = self[split_dy](events, **kwargs)
        events = self[pu_weight](events, **kwargs)
        events = self[muon_weight](events, **kwargs)
        events = self[tau_weight](events, **kwargs) 
           
    # features
    events = self[dilepton_mass](events, **kwargs)
    events = self[mT](events, **kwargs)
    return events