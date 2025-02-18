function result = mm1crossentropy(zlims, origK, origq, L, mu, lambda)

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
    
    lp = P^L;

    cross_entropy = 0;
    for z=[zlims(1):zlims(2)]
        forecast = 0;
        marginal = 0;
        for y=[1:K]
            forecast = forecast + lp(q,y)*conditionals{y}(z);
            marginal = marginal + pim(y)*conditionals{y}(z);
        end
        % Avoid log(0) and handle cases where forecast or marginal are zero.
        if forecast > 0 && marginal > 0
            cross_entropy = cross_entropy - marginal * log(forecast/marginal);
        elseif marginal > 0 && forecast == 0 % forecast is zero, but marginal isn't.  This is maximally "bad"
            cross_entropy = Inf; % or some other large value, since -log(0) is undefined.
            break; %Exit loop, cross entropy is already infinity.
        end
        % If marginal = 0, then the term in the summation is zero (0 * -inf = 0 by L'Hopital).

    end
    result = cross_entropy;
end