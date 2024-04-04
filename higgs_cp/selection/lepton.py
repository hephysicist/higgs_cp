# coding: utf-8

"""
Lepton selection methods.
"""

from __future__ import annotations

from columnflow.selection import Selector, SelectionResult, selector
from columnflow.columnar_util import set_ak_column
from columnflow.util import DotDict, maybe_import
from columnflow.production.util import attach_coffea_behavior



np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")



@selector(
    uses={
        # Muon nano columns
        f"Muon.{var}" for var in [
        "pt", "eta","phi", "dxy", "dz", "mediumId", "pfRelIso04_all",
        ] 
    },
    exposed=False,
)
def study_muon_selection(
    self: Selector,
    events: ak.Array,
    **kwargs,
) -> tuple[ak.Array, ak.Array, ak.Array]:
    
    # Muon selection returning set of selected events and a mask of selected objects.
 
    selections = {
        "muon_pt_26"          : events.Muon.pt > 26,
        "muon_eta_2p4"        : abs(events.Muon.eta) < 2.4,
        "mediumID"            : events.Muon.mediumId == 1,
        "muon_dxy_0p045"      : abs(events.Muon.dxy) < 0.045,
        "muon_dz_0p2"         : abs(events.Muon.dz) < 0.2,
        #"muon_iso_0p15"      : events.Muon.pfRelIso04_all < 0.15 #This variable is used in categories' definition to estimate QCD
    }
    object_mask = ak.ones_like(events.Muon.pt, dtype=np.bool_)
    selection_steps = {}
    for cut_name in selections.keys():
        object_mask = object_mask & selections[cut_name]
        selection_steps[cut_name] = np.array(ak.sum(selections[cut_name], axis=1) > 0, dtype=np.bool_)
    return events, SelectionResult(
        steps = selection_steps,
    ), object_mask
    
@selector(
    uses={
        # Tau nano columns
        f"Tau.{var}" for var in [
        "pt","eta","dz", 
        "idDeepTau2018v2p5VSjet", "idDeepTau2018v2p5VSe", "idDeepTau2018v2p5VSmu", 
        "decayMode", "decayModePNet",
        ] 
    },
    exposed=False,
)

def study_tau_selection(
    self: Selector,
    events: ak.Array,
    **kwargs,
) -> tuple[ak.Array, ak.Array, ak.Array]:
    # """
    # Tau selection returning two sets of indidces for default and veto muons.
   
    # pt sorted indices for converting masks to indices
    deep_tau = self.config_inst.x.deep_tau
    selections = {
        "DeepTauVSjet"  : events.Tau.idDeepTau2018v2p5VSjet >= deep_tau.vs_e_jet_wps[deep_tau.vs_jet], #Medium DeepTau VS jet working point
        "DeepTauVSe"    : events.Tau.idDeepTau2018v2p5VSe   >= deep_tau.vs_e_jet_wps[deep_tau.vs_e], #VVLoose DeepTau VS electron working point
        "DeepTauVSmu"   : events.Tau.idDeepTau2018v2p5VSmu >= deep_tau.vs_mu_wps[deep_tau.vs_mu], #Tight  DeepTau VS muon working point
        "tau_eta_2p3"   : abs(events.Tau.eta) < 2.3,
        "tau_dz_0p2"    : abs(events.Tau.dz) < 0.2,
        "tau_pt_20"     : events.Tau.pt > 20,
        "tau_1and3prong": (events.Tau.decayMode != 5) & (events.Tau.decayMode != 6) #Remove 2 prong decays
        }
    object_mask = ak.ones_like(events.Tau.pt, dtype=np.bool_)
    selection_steps = {}
    #from IPython import embed; embed()
    for cut_name in selections.keys():
        object_mask = object_mask & selections[cut_name]
        selection_steps[cut_name] = np.array(ak.sum(selections[cut_name], axis=1) > 0, dtype=np.bool_)
    return events, SelectionResult(
        steps = selection_steps,
    ), object_mask
    
    
def select_from_multiple_pairs(events: ak.Array,
                        preselection: ak.Array,
                        muons: ak.Array,
                        taus: ak.Array,
                        muon_idx: ak.Array,
                        tau_idx: ak.Array,)-> tuple[selected_events  : ak.Array,
                                                    selected_muon_idx: ak.Array,
                                                     selected_tau_idx: ak.Array]:
                            
    empty_indices = ak.zeros_like(1 * events.event, dtype=np.int32)[:, np.newaxis][..., :0]
    selected_muon_idx = empty_indices
    selected_tau_idx  = empty_indices
    selected_events   = ak.zeros_like(1 * events.event, dtype=np.bool_)
    
    n_multipair_events = ak.sum(ak.sum(preselection, axis=1) > 1)
    if n_multipair_events > 0: 
        #print("Check the multiple pair selection:")                
        muons = ak.drop_none(ak.mask(muons, preselection), behavior=coffea.nanoevents.methods.vector.behavior)
        muons = ak.with_name(muons, "PtEtaPhiMLorentzVector") #Keep Lorentz vector arithmetics
        muon_idx = ak.drop_none(ak.mask(muon_idx, preselection))
        taus = ak.drop_none(ak.mask(taus, preselection), behavior=coffea.nanoevents.methods.vector.behavior)
        taus = ak.with_name(taus, "PtEtaPhiMLorentzVector") #Keep Lorentz vector arithmetics
        tau_idx = ak.drop_none(ak.mask(tau_idx, preselection))
        
        vars4selection = ['muon_pfRelIso04_all', 'muon_pt', 'tau_rawDeepTau2018v2p5VSjet', 'tau_pt']
        is_equal = ak.ones_like(1 * events.event, dtype=np.bool_) #Init is_equal to all Trues allows proper selection during the first loop iteration
        
        for var in vars4selection:
            leptons = muons if "muon" in var else taus
            varname = var.split('_',1)[1]
            do_ascending = (var == 'muon_pfRelIso04_all')# Muons with low iso are the best!
            pair_idx_sorted = ak.argsort(leptons[varname],
                                        axis=1,
                                        ascending=do_ascending) 

            sorted_array = ak.pad_none(leptons[pair_idx_sorted],2, axis=1)
            leading_particle      = sorted_array[:,0]
            subleading_particle   = sorted_array[:,1]
        
            first_muon = muon_idx[pair_idx_sorted][:,:1]
            first_tau  = tau_idx[pair_idx_sorted][:,:1]
            #print(f"{varname}: lead = {leading_particle[ak.num(muons,axis=1)>1][varname][0]}, sublead = {subleading_particle[ak.num(muons,axis=1)>1][varname][0]}")
            if do_ascending:
                selection = ak.fill_none(leading_particle[varname] < subleading_particle[varname],False) #Lower isolation for muon
            else:
                selection = ak.fill_none(leading_particle[varname] > subleading_particle[varname],False) #mu/tau pt and DeepTau score
                
            selection = selection & is_equal # Choose only events that had equal values of the seleciton variable
           
            #print(f"selection: {ak.any(selection)}, is equal: {ak.any(is_equal)}")
            buf_muon_idx = empty_indices
            buf_tau_idx  = empty_indices
            buf_muon_idx = ak.where(selection,first_muon,buf_muon_idx)
            buf_tau_idx  = ak.where(selection,first_tau, buf_tau_idx)
            
            selected_muon_idx = ak.concatenate([buf_muon_idx,
                                                selected_muon_idx],
                                            axis=1)
            selected_tau_idx  = ak.concatenate([buf_tau_idx,
                                                selected_tau_idx],
                                            axis=1)
            selected_events = (selected_events | selection)
            is_equal = (ak.fill_none(leading_particle[varname] == subleading_particle[varname], False) & ~ak.is_none(leading_particle[varname]))
        
    return selected_events, selected_muon_idx, selected_tau_idx
    
                    
@selector(
    uses = 
    {
        f"Muon.{var}" for var in ["pt", "eta","phi", "charge", "dxy", "dz", 
                                  "mediumId", "pfRelIso04_all",]
    } | {
        f"Tau.{var}" for var in ["pt","eta","phi", "dz", "dxy", "charge","rawDeepTau2018v2p5VSjet",
                                 "idDeepTau2018v2p5VSjet", "idDeepTau2018v2p5VSe", "idDeepTau2018v2p5VSmu",] 
    } | {
        f"MET.{var}" for var in ["pt", "phi"]
    } | {
        attach_coffea_behavior
    },
    exposed=False,
)
    

def mutau_selection(
    self: Selector,
    events: ak.Array,
    muon_mask: ak.Array,
    tau_mask : ak.Array,
    **kwargs,
) -> tuple[ak.Array, ak.Array, ak.Array, ak.Array]:  
          
    muon_indices = ak.local_index(events.Muon, axis=1)
    preselected_muons = ak.drop_none(ak.mask(events.Muon, muon_mask),
                                     behavior=coffea.nanoevents.methods.vector.behavior)
    preselected_muons = ak.with_name(preselected_muons, "PtEtaPhiMLorentzVector") #Keep Lorentz vector arithmetics
    preselected_muon_indices = ak.drop_none(ak.mask(muon_indices, muon_mask))
    
    tau_indices = ak.local_index(events.Tau, axis=1)
    preselected_taus = ak.drop_none(ak.mask(events.Tau, tau_mask),
                                    behavior=coffea.nanoevents.methods.vector.behavior)
    preselected_taus = ak.with_name(preselected_taus, "PtEtaPhiMLorentzVector")
    preselected_tau_indices = ak.drop_none(ak.mask(tau_indices, tau_mask)) 
    empty_events = ak.zeros_like(1 * events.event, dtype=np.uint16)
    empty_indices = empty_events[:, np.newaxis][..., :0]
    
    #Produce pairs the preselected muons and taus
    pair_mu, pair_tau = ak.unzip(ak.cartesian([preselected_muons,
                                               preselected_taus], axis=1))
    
    pair_mu_raw_idx, pair_tau_raw_idx = ak.unzip(ak.cartesian([preselected_muon_indices,
                                                     preselected_tau_indices], axis=1))
    #is_os               = (pair_mu.charge * pair_tau.charge) < 0 #Opposite sign mask, I created a custom producer to save this variable and create categories for QCD estimation
    deltaR_selection    = pair_mu.delta_r(pair_tau) > 0.5 #dR between pair constituents
    mT_selection        = pair_mu.mT < 50 # effective mass of muon and MET < 50 GeV/c2
    mutau_selections = {}
    #mutau_selections['mutau_os']        = np.array(ak.sum(is_os,            axis=1) > 0, dtype=np.bool_) #any of mutau pairs of the event has opposite charge
    #mutau_selections['mutau_dr_0p5']    = np.array(ak.sum(deltaR_selection, axis=1) > 0, dtype=np.bool_)
    #mutau_selections['mutau_mt_50']     = np.array(ak.sum(mT_selection,     axis=1) > 0, dtype=np.bool_)
    
    pair_preselection = deltaR_selection & mT_selection #Removed is_os variable to create QCD categories
    
    #Selection of a events with a single pair passed the preselection
    single_pair = (ak.num(pair_preselection, axis=1) == 1) & (ak.sum(pair_preselection, axis=1) == 1) #Check that event contains exactly one pair and this pair passed the preselection
    
    
    
    #Make object-level selections for muons and taus for single pair case
    single_muon_indices = empty_indices
    single_muon_indices = ak.where(single_pair,
                                   pair_mu_raw_idx,
                                   single_muon_indices) 
    single_muon_indices = ak.values_astype(single_muon_indices, np.int32)
    
    single_tau_indices = empty_indices
    single_tau_indices = ak.where(single_pair,
                                  pair_tau_raw_idx,
                                  single_tau_indices) 
    single_tau_indices = ak.values_astype(single_tau_indices, np.int32)
    mutau_selections['single_pair'] = np.array(single_pair , dtype=np.bool_)
    # sel_multiple_pairs, sel_muon_idx, sel_tau_idx = select_from_multiple_pairs(events,
    #                                                         pair_preselection,
    #                                                         pair_mu,
    #                                                         pair_tau,
    #                                                         pair_mu_raw_idx,
    #                                                         pair_tau_raw_idx)
    
    # print(f'selected indices: muons {sel_muon_idx[sel_multiple_pairs]}, taus - {sel_tau_idx[sel_multiple_pairs]}')
    # print(f'N_pres_mu: {ak.num(preselected_muons[sel_multiple_pairs],axis=1)}, \
    #     N_pres_tau: {ak.num(preselected_taus[sel_multiple_pairs],axis=1)}')
    
    # mutau_selections['sel_pairs'] = np.array(sel_multiple_pairs, dtype=np.bool_)
    # if ak.any(sel_multiple_pairs):
    #    pass
    # sel_muon_idx = ak.drop_none(ak.firsts(sel_muon_idx, axis=1))
    # sel_tau_idx = ak.drop_none(ak.firsts(sel_tau_idx, axis=1))
    
    return events, SelectionResult(
        steps=mutau_selections,
    ), single_muon_indices, single_tau_indices,
    
@selector(
    uses={
        # Muon nano columns
        f"Muon.{var}" for var in [
        "pt", "eta","phi", "dz", "dxy", "mediumId", "pfRelIso04_all",
        ] 
    },
    exposed=False,
)
def extra_lepton_veto(
    self: Selector,
    events: ak.Array,
    pair_mu_idx,
    pair_tau_idx,
    **kwargs,
) -> ak.Array :
    
    empty_indices = ak.zeros_like(1 * events.event, dtype=np.int32)[:, np.newaxis][..., :0]
    #Prepare the mask for pair muons
    raw_muon_idx = ak.local_index(events.Muon, axis=1)
    
    (buf_idx, _) = ak.broadcast_arrays(ak.firsts(pair_mu_idx,axis=1)[:,np.newaxis],raw_muon_idx)
    pair_mu_mask = ak.fill_none(buf_idx == raw_muon_idx, False)
    #Select extra muons i.e. those that are not in the pairs
    extra_muons = ak.with_name(ak.drop_none(ak.mask(events.Muon, ~pair_mu_mask),
                                            behavior=coffea.nanoevents.methods.vector.behavior),
                               "PtEtaPhiMLorentzVector")
    #Select pair muons
    pair_muons =  ak.with_name(ak.drop_none(ak.mask(events.Muon, pair_mu_mask),
                                            behavior=coffea.nanoevents.methods.vector.behavior),
                               "PtEtaPhiMLorentzVector")
                               
    #Prepare the mask for pair taus
    raw_tau_idx = ak.local_index(events.Tau, axis=1)
    (buf_idx, _) = ak.broadcast_arrays(ak.firsts(pair_tau_idx,axis=1)[:,np.newaxis],raw_tau_idx)
    pair_tau_mask = ak.fill_none(buf_idx == raw_tau_idx, False)
    #Select pair taus
    pair_taus = ak.with_name(ak.drop_none(ak.mask(events.Tau, pair_tau_mask), 
                                          behavior=coffea.nanoevents.methods.vector.behavior),
                             "PtEtaPhiMLorentzVector")
    #Calculate dR between all extra leptons and pair muons
    (broadcasted_pair_mu, _) = ak.broadcast_arrays(ak.firsts(pair_muons,axis=1)[:,np.newaxis],extra_muons)
    mu_extralep_dr  = ak.fill_none(extra_muons.delta_r(broadcasted_pair_mu), -1)
    #Calculate dR between all exrea leptons and pair taus
    (broadcasted_pair_tau, _) = ak.broadcast_arrays(ak.firsts(pair_taus,axis=1)[:,np.newaxis],extra_muons)
    tau_extralep_dr  = ak.fill_none(extra_muons.delta_r(broadcasted_pair_tau),-1)

    #Perform the selection of extra muons
    extra_lepton_mask = extra_muons.pt > 10
    extra_lepton_mask = extra_lepton_mask & (abs(extra_muons.eta) < 2.5)
    extra_lepton_mask = extra_lepton_mask & (extra_muons.mediumId == 1)
    extra_lepton_mask = extra_lepton_mask &  (extra_muons.pfRelIso04_all < 0.3)
    extra_lepton_mask = extra_lepton_mask & (mu_extralep_dr > 0.5)
    extra_lepton_mask = extra_lepton_mask & (tau_extralep_dr > 0.5)
    
    
    extra_lepton_veto = np.array(~ak.any(extra_lepton_mask, axis=1),dtype=np.bool_) 
    
    pair_mu_idx = ak.where(extra_lepton_veto,
                           pair_mu_idx,
                           empty_indices)
    pair_tau_idx = ak.where(extra_lepton_veto,
                            pair_tau_idx,
                            empty_indices)
    
    return events, SelectionResult(
        steps = {
            "extra_lep_veto" : extra_lepton_veto,
        },
        # objects={
        #     "Muon": {   
        #         "Muon": pair_mu_idx
        #     },
        #     "Tau": {
        #         "Tau": pair_tau_idx
        #     }
        # }
    ), pair_mu_idx, pair_tau_idx
    
    

@selector(
    uses={
        # Muon nano columns
        f"Muon.{var}" for var in [
        "pt", "eta", "phi", "mass", "charge", "mediumId", "pfRelIso04_all",
        "isPFcand", "isGlobal"] #, "isTracker"
    # } | {
    #     f"Tau.{var}" for var in ["pt", "eta", "phi", "dz", "dxy", "charge"] 
    # } | {
    #      f"Electron.{var}" for var in ["pt", "eta", "phi", "dz", "dxy", 
    #                                    "charge", "pfRelIso03_all"] 
    },
    exposed=False,
)
def dilepton_veto(
    self: Selector,
    events: ak.Array,
    pair_mu_idx,
    pair_tau_idx,
    **kwargs,
) -> ak.Array:
    
    empty_indices = ak.zeros_like(1 * events.event, dtype=np.int32)[:, np.newaxis][..., :0]
    
    raw_muon_idx = ak.local_index(events.Muon, axis=1)
    pair_mu1, pair_mu2 = ak.unzip(ak.combinations(events.Muon, 2, axis=1))
    dR_mu1_mu2 = pair_mu1.delta_r(pair_mu2)
    dimu_mask = ak.zeros_like(ak.local_index(pair_mu1, axis=1), dtype=np.bool_)
    
    for muon in [pair_mu1, pair_mu2]:
        dimu_mask = dimu_mask & (muon.pt > 15)
        dimu_mask = dimu_mask & (abs(muon.eta) < 2.4)
        dimu_mask = dimu_mask & muon.isGlobal
        #dimu_mask = dimu_mask & muon.isTracker
        dimu_mask = dimu_mask & muon.isPFcand
        dimu_mask = dimu_mask & (muon.dz < 0.2)
        dimu_mask = dimu_mask & (muon.dxy < 0.045)
        dimu_mask = dimu_mask & (muon.pfRelIso04_all < 0.3)
    
    dimu_mask = dimu_mask & ((pair_mu1.charge * pair_mu2.charge) < 0)
    dimu_mask = dimu_mask & (dR_mu1_mu2 > 0.15)
    dilepton_veto = np.array(~ak.any(dimu_mask, axis=1),dtype=np.bool_) 
    #print(ak.sum(ak.num(raw_muon_idx)), ak.sum(ak.sum(dimu_mask,axis=1)))
    return events, SelectionResult(
        steps = {
            "dilep_veto" : dilepton_veto,
        }
    )
