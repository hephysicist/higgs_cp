import functools
from columnflow.production import Producer, producer
from columnflow.util import maybe_import
from columnflow.columnar_util import set_ak_column, EMPTY_FLOAT, Route
from columnflow.production.util import attach_coffea_behavior
from IPython import embed
ak = maybe_import("awkward")
np = maybe_import("numpy")

# helper
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)


@producer(
    uses = 
    {
        f"Muon.{var}" for var in ["pt", "eta","phi", "mass","charge"]
    } | {
        f"Tau.{var}" for var in ["pt","eta","phi", "mass", "dxy", "dz", "charge"] 
    } | {attach_coffea_behavior},
    produces={
        "mutau_mass"
    },
)
def dilepton_mass(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    print("Producing dilepton mass...")
    events = self[attach_coffea_behavior](events, **kwargs)
    
    muon = ak.firsts(events.Muon, axis=1)
    tau = ak.firsts(events.Tau, axis=1)
    mutau_obj = muon + tau
    mutau_mass = ak.where(mutau_obj.mass2 >=0, mutau_obj.mass, EMPTY_FLOAT)
    events = set_ak_column_f32(events,"mutau_mass",mutau_mass)
    
    return events

@producer(
    uses = 
    {
        f"Muon.{var}" for var in ["pt","phi"]
    } | {
        f"MET.{var}" for var in ["pt","phi"] 
    } | {attach_coffea_behavior},
    produces={
        "Muon.mT"
    },
)
def mT(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    
    # def calc_mt(muon, met):
    #     cos_dphi = np.cos(muon.phi - met.phi)
    #     mT_vals = np.sqrt(2*muon.pt*met.pt*(1. - cos_dphi))
    #     return ak.fill_none(mT_vals, EMPTY_FLOAT)
    
    # initial_dim = ak.num(events.Muon)
    # (_, met_br) = ak.broadcast_arrays(events.Muon, events.MET[:,np.newaxis])
   
    # flat_muon = ak.flatten(events.Muon, axis=1)
    # flat_met  = ak.flatten(met_br,      axis=1)
    #  mT_values = ak.unflatten(calc_mt(flat_muon, flat_met), initial_dim,axis=0)
    cos_dphi = np.cos(events.Muon.delta_phi(events.MET))
    mT_values = np.sqrt(2 * events.Muon.pt * events.MET.pt * (1 - cos_dphi))
    mT_values = ak.fill_none(mT_values, EMPTY_FLOAT)
    events = set_ak_column_f32(events, Route("Muon.mT"), mT_values)
    return events
    

    
   