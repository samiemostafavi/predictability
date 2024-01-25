import numpy as np

def transition_matrix_power(transition_matrix, power):
    return np.linalg.matrix_power(transition_matrix, power)

class MarkovChain:
    def __init__(self,transition_matrix):
        self.transition_matrix = transition_matrix
        self.K = len(transition_matrix[0])

    def transition_matrix_power(self, power):
        return np.linalg.matrix_power(self.transition_matrix, power)

    def get_stationary_probs(self, power = 1000):
        initial_state = np.zeros(len(self.transition_matrix[0]))
        initial_state[0] = 1
        new_r = transition_matrix_power(self.transition_matrix,power)
        res = initial_state.dot(new_r)
        return res
    
    def simulate(self, initial_state, num_steps):
        num_states = len(self.transition_matrix)
        states = [initial_state]

        for _ in range(num_steps - 1):
            current_state = states[-1]
            transition_probs = self.transition_matrix[current_state, :]
            next_state = np.random.choice(np.arange(num_states), p=transition_probs)
            states.append(next_state)

        return states