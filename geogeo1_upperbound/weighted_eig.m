function [sortedEigenvectors, sortedEigenvalues] = weighted_eig(P, pim)
    [m, n] = size(P);
    if m ~= n
        error('Matrix P must be square');
    end
    if length(pim) ~= n
        error('Vector pim must be the same size as matrix P dimensions');
    end
    

    B = diag(sqrt(pim));
    Binv = diag(sqrt(pim.^-1));
    
    A = B*P*Binv;
%     if ~issymmetric(A)
%         A = Binv*P*B;
%         if ~issymmetric(A)
%             error('A is not symmetric, probably P is not reversible');
%         end
%     end
    
    [V, D] = eig(A);

    eigenvalues = diag(D);
    eigenvectors = Binv*V;
    
    % Sort eigenvalues in ascending order
    [sortedEigenvalues, idx] = sort(eigenvalues,'descend');

    % Reorder eigenvectors accordingly
    sortedEigenvectors = eigenvectors(:, idx);
    
    % Make sure column 1 in eigenvectors is (1,1,...,1) not (-1,-1,...,-1)
    if sortedEigenvectors(1,1) < 0
        sortedEigenvectors = -1*sortedEigenvectors;
    end
    
    sortedEigenvalues = diag(sortedEigenvalues);
end
