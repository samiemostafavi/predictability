function result = predictability_exact_multi(zlims, qs, L, Ps, pims, M, Ks, conditionals)
    forecast_conv = [];
    marginal_conv = [];
    forecast_dists = [];
    marginal_dists = [];
    for idx=[1:M]
        P = Ps{idx};
        q = qs(idx);
        K = Ks(idx);
        pim = pims(idx,:);
        lp = P^L;
        forecast_dist = [];
        marginal_dist = [];
        for z=[zlims(1):zlims(2)]
            forecast = 0;
            marginal = 0;
            for y=[1:K]
                marginal = marginal + pim(y)*conditionals{y,idx}(z); % must start from zero
                if q >= 0
                    forecast = forecast + lp(q,y)*conditionals{y,idx}(z); % must start from zero
                else % if q is negative it means we have no observation
                    forecast = marginal;
                end
                
            end
            forecast_dist = [forecast_dist, forecast];
            marginal_dist = [marginal_dist, marginal];
        end
        forecast_dists = [forecast_dists; forecast_dist];
        marginal_dists = [marginal_dists; marginal_dist];
        if ~isempty(forecast_conv)
            forecast_conv = conv(forecast_conv,forecast_dist);
        else
            forecast_conv = forecast_dist;
        end
        if ~isempty(marginal_conv)
            marginal_conv = conv(marginal_conv,marginal_dist);
        else
            marginal_conv = marginal_dist;
        end  
    end
    result = 1/2*sum(abs(forecast_conv-marginal_conv));
end
