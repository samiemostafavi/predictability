import numpy as np
from scipy.stats import norm

from predictability.gmm import GaussianMixtureModel,GMMComponent
from predictability.mc import MarkovChain
from predictability.core import calc_forecast_weights

if __name__ == '__main__':
    mc = MarkovChain(
        np.array([
            [0.9, 0.05, 0.05],
            [0.25, 0.5, 0.25],
            [0.15, 0.25, 0.6]
        ])
    )
    
    # define prior distributions
    # analysis
    analysis_dists = []
    for k in range(mc.K):
        analysis_dists.append(norm(loc=(k+1)*3, scale=1))

    # estimate the climatological dist. which is a mixture
    clim_weights = mc.get_stationary_probs()
    clim_gmm = GaussianMixtureModel(
        weights = clim_weights, 
        components=[
            GMMComponent(dist.mean(), dist.std()) for dist in analysis_dists
        ]
    )
