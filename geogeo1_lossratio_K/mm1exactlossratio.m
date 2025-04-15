function result = mm1exactlossratio(origK, origq, L, mu, lambda)

    K = origK;
    q = origq;
    % form P matrix
    P = sparse(K, K);
    P(1, 1) = 1 - lambda * (1 - mu);
    P(sub2ind([K, K], 2:K, 1:K-1)) = lambda * (1 - mu);
    P(sub2ind([K, K], 1:K-1, 2:K)) = mu * (1 - lambda);
    P(sub2ind([K, K], 2:K-1, 2:K-1)) = lambda * mu + (1 - lambda) * (1 - mu);
    P(K, K) = 1 - mu * (1 - lambda);
    P = full(P)';
    
    tmp = P^10000;
    pim = tmp(1,:);
    
    lp = P^L;
    forecast = lp(q,K);
    marginal = pim(K);
    result = abs(forecast-marginal);
end
