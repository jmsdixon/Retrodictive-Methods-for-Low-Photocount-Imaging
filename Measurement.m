% Author: James Dixon
% Simulated measurement program.  Returns an effective frame of any size,
% 'eff_f_sz', from a simulated Monte Carlo simulated measurement by a SPAD 
% array of a coherent source given an input image, 'n_bars_in'.  The input
% is a matrix of mean photon numbers incident on each detector.
% 'sz1' is number of rows of SPADs
% 'sz2' number of columns of SPADs
% 'fpf' number of frames-per-frame. set to 2 for reported work

function [measurements,eff_frame]=Measurement(n_bars_in,eff_f_sz,fpf,sz1,sz2)
eta = 0.5;
measurements = zeros([sz1 sz2 eff_f_sz]);

dist = @(nbar)exp(-eta*nbar); % No click probabilities, coherent source
nc_probs = dist(n_bars_in);   % Matrix - no click prob at each detector

for i=1:1:eff_f_sz
    % Generate photo-counts if rand num >= no click probability
    measurements(:,:,i) = sum((rand([sz1 sz2 fpf]) >= nc_probs),3);
end

eff_frame = sum(measurements,3);
end