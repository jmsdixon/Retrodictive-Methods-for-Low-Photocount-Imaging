% Author: James Dixon
% Training 2 variable PMF for n_bars
close all; clear;
%% 
% User set specifications for training n_bar probability mass function PMF 
% P(Mj, Mjo | n_bar)
%
% Generates 'n_bar_dists' with all Mj and Mjo counts for each discrete 
% n_bar and normalises them to get 'dists_norm' PMF
%
% Set parameters and training image, then run for say 50 sets. Check
% 'timer' before running longer.
% Check distributions are collecting all counts,
% i.e. that 'max_Mj_counts' is large enough.  If some n_bars are missing, 
% as indicated in 'missing_n_bars', training on several images is the
% best way to ensure all n_bars are represented.

%% User sets paramters here

frames_per_frame = 2; % Set to 2 for Dr Johnstone's lab data
eta = 0.5;

% Set n_bar_max
n_bar_max = 0.0255;    % 0.0255 works well for with Johnstone's lab data
lstep = n_bar_max/255;
round_to = 4;

sz1 = 168; % image height in pixels. Set to correspond to 'Test-target' size
sz2 = 128; % image width in pixels

max_Mj_counts = 8; % Size of each square pmf per n_bar, including zero.
                   % 8 works well for individual frames but needs to be
                   % bigger for effective frames of size > 1 frame. Best to
                   % run for 50 sets and see if distributions are not
                   % being limited. If too small, leads to algorithm
                   % failure when imaging.

eff_f_sz = 1;     % Effective frame size for imaging algorithm
data_sets = 10;   % Number of data sets for statistics collection. 
                  % A minimum of 5000 sets for the full PMF is recommended


dist = @(nbar)exp(-eta*nbar); % No click probabilities, coherent source

%% User set training image PNG file
%   Open and format training image as 'n_bars_train'.

fname='test_target.png'; % Set png training image
train_im=imread(fname); train_im=rgb2gray(train_im); 

train_im=imresize(train_im,[sz1 sz2]);
n_bars_train=double(train_im);

n_bars_train=n_bars_train.*(n_bar_max/255);
n_bars_train=round(n_bars_train,round_to);

%% Some storage arrays

n_bars = 0:lstep:n_bar_max;
n_bars = round(n_bars,round_to);
n_bar_dists = zeros(max_Mj_counts,max_Mj_counts,length(n_bars)); 
% 'n_bar_dists' is non-normalised PMF

%% Finding which n_bars are missing if any

missing_n_bars = zeros(1,length(n_bars));
idx = 1;
for i=0:lstep:n_bar_max
    is_in = round(i,round_to)==n_bars_train;
    if sum(is_in,'all') == 0
        missing_n_bars(idx) = 1;
    end
    idx = idx + 1;
end

keepers = missing_n_bars==0;
missing_n_bars = n_bars.*missing_n_bars;
n_bar_keepers = keepers.*n_bars;

%% Use input mean n_bars to generate counts and store Mj and Mjo per n_bar

in_nc_probs = dist(n_bars_train); % Matrix - no click prob at each detector
tic;
for data_set=1:1:data_sets
    % If rand nums are >= no click probability, then a count is produced
    measurements = sum((rand([[sz1 sz2] frames_per_frame*eff_f_sz]) >= in_nc_probs),3);


    %% Get smoothing algorithmed matrices

    RS = ones(sz1,sz2);
    R0 = 1;
    [Mj,Nj,Mj0,Nj0,Mbar]=Bayret_mat_Mj(measurements,RS,R0);


    %% Find Mj and MJ0 statistics for each keeper n_bar
    Mj_plus1 = Mj + 1;
    Mj0_plus1 = Mj0 + 1;

    for nb_idx=1:1:length(n_bars)
    
        if nb_idx==1 % For n_bar=0
            targets = n_bars_train == 0;  % 
            for Mj_count=1:1:max_Mj_counts

                    targMj = targets.*Mj_plus1 == Mj_count; % Mj elements where this n_bar got counts
                    targMj0 = targMj.*Mj0_plus1;           % Mj0 elements where this n_bar Mj got counts

                for Mj0_count=1:1:max_Mj_counts

                    n_bar_dists(Mj_count,Mj0_count,nb_idx) = n_bar_dists(Mj_count,Mj0_count,nb_idx) + sum((targMj0==Mj0_count),'all');

                end
            end
        end
        
        if n_bar_keepers(nb_idx)==0  % Where n_bar is not in n_bars_train
            continue                 %  and zero
        else
            n_bar = n_bar_keepers(nb_idx);
        end

        targets = n_bars_train == round(n_bar,round_to);  % round n_bar to deal with NaNs
        for Mj_count=1:1:max_Mj_counts

                targMj = targets.*Mj_plus1 == Mj_count; % Mj elements where this n_bar got counts
                targMj0 = targMj.*Mj0_plus1;           % Mj0 elements where this n_bar Mj got counts

            for Mj0_count=1:1:max_Mj_counts

                n_bar_dists(Mj_count,Mj0_count,nb_idx) = n_bar_dists(Mj_count,Mj0_count,nb_idx) + sum((targMj0==Mj0_count),'all');

            end
        end
    end
end
timer=toc;

%% Normalise n_bar Mj count distributions
dists_norm = bsxfun(@rdivide,n_bar_dists,sum(n_bar_dists,[1 2]));


%% Setting zero elements to 1e-15
% This prevents Nans appearing in retrodictions
inds = dists_norm==0;
dists_norm(inds) = 1e-15;
