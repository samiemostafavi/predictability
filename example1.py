import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt

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

    # plot analysis and climatological distributions
    # analysis
    samples_per_dist = 10000
    plt.figure(figsize=(6, 3))
    samples = [dist.rvs(size=samples_per_dist) for dist in analysis_dists]
    for k in range(mc.K):
        plt.hist(samples[k], bins=100, density=True, alpha=0.5, label=f'State {k}')

    plt.xlabel('Value')
    plt.ylabel('Probability Density')
    plt.legend()
    plt.savefig("res1.jpg")

    # climatological
    num_samples = 10000
    plt.figure(figsize=(6, 3))
    samples = clim_gmm.sample(num_samples)
    plt.hist(samples, bins=100, density=True, alpha=0.5, label=f'Climatalogical')

    plt.xlabel('Value')
    plt.ylabel('Probability Density')
    plt.legend()
    plt.savefig("res2.jpg")


    Ts = np.array(list(range(40)))+1
    initial_state = 0

    # climatological entropy calc
    clim_lb_entropy = clim_gmm.calc_entropy_lb()
    clim_approx_entropy = clim_gmm.calc_entropy_approx()
    clim_ub_entropy = clim_gmm.calc_entropy_ub()

    # forecast entropy calc
    forecast_lb_entropies = []
    forecast_approx_entropies = []
    forecast_ub_entropies = []
    for T in Ts:
        new_weights = calc_forecast_weights(mc, initial_state, T)
        forecast_gmm = GaussianMixtureModel(
            weights = new_weights, 
            components=[
                GMMComponent(dist.mean(), dist.std()) for dist in analysis_dists
            ]
        )
        forecast_lb_entropies.append(forecast_gmm.calc_entropy_lb())
        forecast_approx_entropies.append(forecast_gmm.calc_entropy_approx())
        forecast_ub_entropies.append(forecast_gmm.calc_entropy_ub())

    # Plot entropy
    plt.figure(figsize=(6, 3))
    plt.plot(Ts,clim_lb_entropy*np.ones(len(forecast_lb_entropies)),marker=".",label="Climatological LB Entropy")
    plt.plot(Ts,forecast_lb_entropies,marker=".",label="Forecast LB Entropy")

    plt.plot(Ts,clim_approx_entropy*np.ones(len(forecast_approx_entropies)),marker=".",label="Climatological Approx Entropy")
    plt.plot(Ts,forecast_approx_entropies,marker=".",label="Forecast Approx Entropy")

    plt.plot(Ts,clim_ub_entropy*np.ones(len(forecast_ub_entropies)),marker=".",label="Climatological UB Entropy")
    plt.plot(Ts,forecast_ub_entropies,marker=".",label="Forecast UB Entropy")

    plt.title('Entropy Bounds')
    plt.xlabel('Value')
    plt.ylabel('Entropy')
    plt.legend()
    plt.savefig("res3.jpg")


    # Plot predictability
    plt.figure(figsize=(6, 3))
    plt.plot(Ts,clim_lb_entropy*np.ones(len(forecast_lb_entropies)) - np.array(forecast_lb_entropies),marker=".",label="Lower Bound")
    plt.plot(Ts,clim_approx_entropy*np.ones(len(forecast_approx_entropies)) - np.array(forecast_approx_entropies),marker=".",label="Approx.")
    plt.plot(Ts,clim_ub_entropy*np.ones(len(forecast_ub_entropies)) - np.array(forecast_ub_entropies),marker=".",label="Upper Bound")
    plt.title('Predictability')
    plt.xlabel('Value')
    plt.ylabel('Predictability')
    plt.legend()
    plt.savefig("res4.jpg")