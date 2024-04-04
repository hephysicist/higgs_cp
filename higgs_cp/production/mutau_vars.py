import functools
from columnflow.production import Producer, producer
from columnflow.util import maybe_import
from columnflow.columnar_util import set_ak_column, EMPTY_FLOAT, Route, optional_column as optional
from columnflow.production.util import attach_coffea_behavior
#from IPython import embed
ak = maybe_import("awkward")
np = maybe_import("numpy")
coffea = maybe_import("coffea")
# helper
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)
# @producer(
#     uses = 
#     {
#         f"Muon.{var}" for var in ["pt", "eta","phi", "mass","charge"]
#     } | {
#         f"Tau.{var}" for var in ["pt","eta","phi", "mass", "dxy", "dz", "charge"] 
#     } | {attach_coffea_behavior},
#     produces={
#         "mutau_mass"
#     },
# )
# def dilepton_mass(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
#     print("Producing dilepton mass...")
#     events = self[attach_coffea_behavior](events, **kwargs)
#     muon = ak.firsts(events.Muon, axis=1)
#     tau = ak.firsts(events.Tau, axis=1)
#     mutau_obj = muon + tau
#     mutau_mass = ak.where(mutau_obj.mass2 >=0, mutau_obj.mass, EMPTY_FLOAT)
#     events = set_ak_column_f32(events,"mutau_mass",mutau_mass)
    
#     return events

@producer(
    uses = 
    {
        f"Muon.{var}" for var in ["pt", "eta","phi", "mass","charge"]
    } | {
        f"Tau.{var}" for var in ["pt","eta","phi", "mass", "dxy", "dz", "charge"] 
    } | {optional("Tau.pt_no_tes"), optional("Tau.mass_no_tes")
    } | {attach_coffea_behavior},
    produces={
        "mutau_mass_no_tes", "mutau_mass"
    },
)
def dilepton_mass(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    print("Producing dilepton mass...")
    from coffea.nanoevents.methods import vector
    events = self[attach_coffea_behavior](events, **kwargs)
    muon = ak.firsts(events.Muon, axis=1)
    if self.dataset_inst.is_mc: tags = ["", "_no_tes"]
    else: tags = [""]
    for tag in tags:
        tau = ak.zip(
            {
                "pt": events.Tau[f"pt{tag}"],
                "eta": events.Tau.eta,
                "phi": events.Tau.phi,
                "mass": events.Tau[f"mass{tag}"],
            },
            with_name="PtEtaPhiMLorentzVector",
            behavior=vector.behavior,
        )
        mutau_obj = muon + tau
        mutau_mass = ak.flatten(ak.where(mutau_obj.mass2 >=0, mutau_obj.mass, EMPTY_FLOAT))
        events = set_ak_column_f32(events,f"mutau_mass{tag}",mutau_mass)
    if self.dataset_inst.is_data:
        events = set_ak_column_f32(events,f"mutau_mass_no_tes",mutau_mass)
    return events


@producer(
    uses = 
    {
        "Muon.charge","Tau.charge"
    } | {attach_coffea_behavior},
    produces={
        "rel_charge"
    },
)
def rel_charge(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    print("Producing  pair relative charge...")
    events = self[attach_coffea_behavior](events, **kwargs)
    lep1_charge = ak.firsts(events.Muon.charge, axis=1)
    lep2_charge = ak.firsts(events.Tau.charge, axis=1)
    events = set_ak_column_f32(events,f"rel_charge", lep1_charge * lep2_charge)
    return events

@producer(
    uses = 
    {
        f"Muon.{var}" for var in ["pt","phi"]
    } | {
        f"PuppiMET.{var}" for var in ["pt","phi"] 
    } | {attach_coffea_behavior},
    produces={
        "Muon.mT"
    },
)
def mT(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    print("producing mT...")
    events = self[attach_coffea_behavior](events, **kwargs)
    cos_dphi = np.cos(events.Muon.delta_phi(events.PuppiMET))
    mT_values = np.sqrt(2 * events.Muon.pt * events.PuppiMET.pt * (1 - cos_dphi))
    mT_values = ak.fill_none(mT_values, EMPTY_FLOAT)
    events = set_ak_column_f32(events, Route("Muon.mT"), mT_values)
    return events
    

    
   