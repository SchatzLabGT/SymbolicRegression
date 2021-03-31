function Xi = find_coeffs(Theta, gamma)
% compute sparse regression on Theta * Xi = 0
% first column is (usually) dx/dt term
% gamma : sparsification parameter - higher gamma means more pruning
    [h,w] = size(Theta);
    [~,~,V] = svd(Theta);
    Xi = V(:,end); % first estimate of coefficients (not sparse)
    lambda = norm(Theta*Xi); % residual as threshold
    smallinds = zeros(w,1);
    for i=1:w
          % product of the coefficient and characteristic size of library function
          for p_ind = 1:w
            product(p_ind) = norm(Xi(p_ind)*Theta(:,p_ind)./sqrt(sum(Theta(:,:).^2,2)));
          end
          product(smallinds==1)=Inf;
          [Y,I] = min(product); % find next term to eliminate        
          smallinds(I) = 1;
          if sum(smallinds==0)==0 % no more terms to remove
              break;
          end
          Xi_old = Xi;
          Xi(smallinds==1) = 0; % set coefficients of negligible terms to 0
          [~,~,V] = svd(Theta(:,smallinds==0));
          Xi(smallinds==0) = V(:,end); % find new coefficients
          lambda_old = lambda;
          lambda = norm(Theta*Xi);
          if (lambda/lambda_old > gamma) || nnz(Xi)==0 % condition to end iteration
             Xi = Xi_old;
             break
          end
    end
end