function [aggP, aggConditionals] = aggregateMMP(P, pim, conditionals, Kp)
    
    % form aggregated state transition
    aggP = aggregateMC(P, pim, Kp);
    
    % form aggregated conditional dists.
    K = length(pim);

    % aggregate to Kp even states
    % example for K=10, is Kp = 5
    aggConditionals = cell(Kp,1);
    for a = [1:Kp]
        oa = [(a-1)*K/Kp+1:a*K/Kp];
        aggConditionals{a} = @(z) aggcondsum(z,conditionals,oa,pim);
    end
end

function res = aggcondsum(z,conditionals,oa,pim)
    sumprob = 0;
    aggcondsum = 0;
    for j_ind = [1:length(oa)]
        j = oa(j_ind);
        sumprob = sumprob + pim(j);
        aggcondsum = aggcondsum + pim(j)*conditionals{j}(z);
    end
    res = aggcondsum/sumprob;
end

function result = aggregateMC(P, pim, Kp)

    K = length(pim);
    % aggregate to Kp even states
    % example for K=10, is Kp = 5
    Pa = zeros(Kp);
    for a = [1:Kp]
        oa = [(a-1)*K/Kp+1:a*K/Kp];
        for b = [1:Kp]
            ob = [(b-1)*K/Kp+1:b*K/Kp];
            if a ~= b
                numerator = 0;
                for j_ind = [1:length(oa)]
                    j = oa(j_ind);
                    for i_ind = [1:length(ob)]
                        i = ob(i_ind);
                        numerator = numerator + pim(j)*P(j,i);
                    end
                end
                denumerator = 0;
                for j_ind = [1:length(oa)]
                    j = oa(j_ind);
                    denumerator = denumerator + pim(j);
                end
                Pa(a,b) = numerator/denumerator;
            else
                numerator = 0;
                for j_ind = [1:length(oa)]
                    j = oa(j_ind);
                    for i_ind = [1:length(oa)]
                        i = oa(i_ind);
                        numerator = numerator + pim(j)*P(j,i);
                    end
                end
                denumerator = 0;
                for j_ind = [1:length(oa)]
                    j = oa(j_ind);
                    denumerator = denumerator + pim(j);
                end
                Pa(a,b) = numerator/denumerator;
            end
        end
    end
    result = Pa;


end


