% Author: James Dixon, 2021
% Bayes Theorem using trained stats - 2 variable pmf
% Takes 'Mj' and 'Mj0' counts, the 'prior' mean photon number distribution 
% and trained distridutions, 'dists'.
% 'lmax' is maximum mean photon number
% 'lstep' is the mean photon number discretisation
% 'retro' is the updated posterior to be returned

function [retro]=STFFR(lstep,lmax,prior,Mj,Mj0,dists)

num_dets = 21504;
len_prior = round(lmax/lstep + 1);
preconds = zeros(num_dets,len_prior);

% Reshape input data
Mj = reshape(Mj,num_dets,[])+1;
Mj0 = reshape(Mj0,num_dets,[])+1;

% Retrieve appropriate predictive conditionals from trained stats, dists
for j=1:1:num_dets
    preconds(j,:) = dists(Mj(j),Mj0(j),:);
end
% Compute posterior, retro, with Bayes' Theorem
prods = prior.*preconds;
retro = bsxfun(@rdivide,prods,sum(prods,2));

end