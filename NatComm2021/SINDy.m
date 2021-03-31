%% --------------------------------------------------------------- SINDy()
function Xi = SINDy (Theta, dXdt)
%{
Compute sparse regression on dX = Theta * Xi
Regression technique used: sequential least squares

Modified procedure from:
S. H. Rudy, S. L. Brunton, J. L. Proctor, J. N. Kutz, Data-driven 
discovery of partial differential equations. Sci. Adv. 3, e1602614 (2017)
%}

Xi = Theta \ dXdt;

gamma = 0.05;
lambda = gamma*mean(abs(dXdt)); 

for i = 1:5

  product = zeros(size(Xi)); 
  [~,w] = size(Theta);
  for p_ind = 1:w
    product(p_ind) = mean(abs(Xi(p_ind)*Theta(:,p_ind)));
  end

  smallinds = product < lambda;
  Xi(smallinds) = 0;                % set negligible terms to 0
  for ind = 1:size(dXdt,2)   
    biginds = ~smallinds(:,ind);
    Xi(biginds,ind) = Theta(:,biginds) \ dXdt(:,ind);
  end
end
    
end