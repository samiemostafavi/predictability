import numpy as np
from scipy.stats import norm
import math

class GaussianMixtureModel:
    def __init__(self, weights, components):
        self.weights = weights
        self.components = components
        self.K = len(self.components)

    def evaluate_gradient(self, X):
        gradient_value = 0
        for k in range(self.K):
            gradient_value += self.weights[k] * self.components[k].evaluate_gradient(X)
        return gradient_value

    def compute_H0(self):
        H0 = 0
        for k in range(self.K):
            H0 -= self.weights[k] * self.logpdf(self.components[k].mu)
        return H0

    def compute_H2(self):
        H2 = 0
        for k in range(self.K):
            s = np.sum(np.sum(self.compute_F(self.components[k].mu) * self.components[k].sigma))
            H2 -= self.weights[k] * 0.5 * s
        return H2

    def compute_F(self, X):
        F = 0.0
        inverse_of_pdf = 1.0 / self.pdf(X)
        gradient_in_X = self.evaluate_gradient(X)

        for j in range(self.K):
            lambda_val = 1.0/self.components[j].sigma
            mean_dif = (X - self.components[j].mu)
            F += self.weights[j] * lambda_val * (
                    inverse_of_pdf * ((mean_dif) * gradient_in_X) + (mean_dif) * (lambda_val * mean_dif) - 1.0
            ) * self.components[j].pdf(X)

        F *= inverse_of_pdf
        return F

    def logpdf(self, x):
        # Implement the logpdf function as needed for your application
        return np.log(self.pdf(x))

    def pdf(self, x):
        # Implement the pdf function as needed for your application
        return self.weights @ np.array([component.pdf(x) for component in self.components])

    def sample(self,num_samples):
        return np.concatenate([
            norm(loc = comp.mu,scale = comp.sigma).rvs(size=int(weight * num_samples)) for comp, weight in zip(self.components, self.weights)
        ])

    def calc_entropy_approx(self):
        return self.compute_H0()+self.compute_H2()

    # lower bound on the entropy of a GMM
    # https://ieeexplore.ieee.org/document/4648062 Theorem 2
    def calc_entropy_lb(self):
        outer_sum = 0
        for i,outer_dist in enumerate(self.components):
            w_i = self.weights[i]
            mu_i = outer_dist.mu
            inner_sum = 0
            for j,inner_dist in enumerate(self.components):
                w_j = self.weights[j]
                mu_j = inner_dist.mu
                z_ij = norm.pdf(mu_i, mu_j, math.sqrt(inner_dist.sigma**2+outer_dist.sigma**2))
                inner_sum += w_j*z_ij

            outer_sum += -w_i*math.log(inner_sum)

        return outer_sum

    # upper bound on the entropy of a GMM
    # https://ieeexplore.ieee.org/document/4648062 equation 8
    def calc_entropy_ub(self):
        sum = 0
        for i,dist in enumerate(self.components):
            w_i = self.weights[i]
            if w_i > 0:
                C_i = dist.sigma**2
                sum += w_i*(-math.log(w_i) + 0.5*math.log((2*math.pi*math.e)*C_i))
        return sum


class GMMComponent:
    def __init__(self, mu, sigma):
        self.mu = mu
        self.sigma = sigma

    def evaluate_gradient(self, X):
        return -(1.0/self.sigma) * (X - self.mu) * self.pdf(X)

    def pdf(self, X):
        # Implement the pdf function for a component as needed
        return norm.pdf(X, self.mu, self.sigma)