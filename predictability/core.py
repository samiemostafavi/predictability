import numpy as np

def calc_forecast_weights(markov_chain, initial_state, tau):
    # calculate forecast dist. loc and scale
    forecast_weigths = []
    for k in range(markov_chain.K):
        new_r = markov_chain.transition_matrix_power(tau)
        sum_prob =new_r[initial_state,k]
        forecast_weigths.append(sum_prob)

    forecast_weigths = np.array(forecast_weigths)
    return forecast_weigths