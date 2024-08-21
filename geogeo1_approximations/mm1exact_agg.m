function result = mm1exact_agg(level, zlims, origK, origq, L, mu, lambda)

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
    
    % mm1 conditionals
    conditionals = cell(K,1);
    for n=[1:K]
        conditionals{n} = @(z) nbinpdf(z,n,mu);
    end
    
    % aggregate states
    [aggP, aggConditionals] = aggregateMMP(P, pim, conditionals, origK/level);
    P = aggP;
    q = origq/level;
    conditionals = aggConditionals;
    K = length(aggP);
    tmp = aggP^10000;
    pim = tmp(1,:);
    lp = P^L;

    l1norm = 0;
    for z=[zlims(1):zlims(2)]
        forecast = 0;
        marginal = 0;
        for y=[1:K]
            forecast = forecast + lp(q,y)*conditionals{y}(z);
            marginal = marginal + pim(y)*conditionals{y}(z);
        end
        l1norm = l1norm + abs(forecast-marginal);
    end
    result = 1/2*l1norm;
end
