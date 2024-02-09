"""
Wrappers for some default sets of producers.
"""

from columnflow.production import Producer, producer
from columnflow.util import maybe_import

from columnflow.production.normalization import normalization_weights
from columnflow.production.categories import category_ids
from higgs_cp.production.example import features, cutflow_features 
from higgs_cp.production.mutau_vars import dilepton_mass, mT
from IPython import embed

ak = maybe_import("awkward")


@producer(
    uses={
        category_ids, features, normalization_weights, cutflow_features, dilepton_mass, mT
    },
    produces={
        category_ids, features, normalization_weights, cutflow_features, dilepton_mass, mT
    },
)
def default(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    # category ids
    #events = self[category_ids](events, **kwargs)
    events = self[normalization_weights](events, **kwargs)
    # features
    #embed()
    events = self[features](events, **kwargs)
    events = self[dilepton_mass](events, **kwargs)
    events = self[mT](events, **kwargs)
    #embed()
    return events