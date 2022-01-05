% Author: James Dixon, 2021
% Based on my own work, "How Many Photons Make an Image?",
% 4th year project style retrodictor; optimised.
% 'Click' or 'no click' detection in 'frm' used to update posterior
% distributions at each pixel.  Assumes coherent light source.
% 'lmax' is maximum mean photon number
% 'lstep' is the mean photon number discretisation
% 'prior' is the prior mean photon number distribution to be updated by
% Bayes' Theorem in this function.
% 'retro' is the posterior to be returned.

% frames per frame, 'fpf', variable set to 2 for test_target data from lab.
% Otherwise ~ 30-40 for good results - not reported.

function [retro]=FFR(lstep,lmax,prior,frm,fpf)

num_dets = 21504;
eta=0.5;
len_prior = lmax/lstep + 1;
n_bars=0:lstep:lmax;

% Coherent State: prob of no click given m mean photon number
nc_prob = @(m)exp(-eta*m); % No click prob

% predictive probabilities for click/ no click given n_bar coherent source
predictives=zeros(2,len_prior);
predictives(1,:)=nc_prob(n_bars);
predictives(2,:) = 1 - predictives(1,:);

% Restructure input frame data frm to act as indices.
frm = reshape(frm,num_dets,[]);
frms = repmat(frm,1,fpf);
frms = bsxfun(@minus,frms,0:1:fpf-1);
inds = frms>0;

% Retrodictivelty update posterior intensity distribution
for i=1:1:fpf
    % Compute posterior, retro, with Bayes' Theorem
    preconds = predictives(inds(:,i)+1,:);
    prods = prior.*preconds;
    retro = bsxfun(@rdivide,prods,sum(prods,2));
    prior = retro;
    

end

end