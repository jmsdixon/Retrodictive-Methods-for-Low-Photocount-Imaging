% Author: James Dixon, 2021
% Data collection and user interface script
close all; clear;

%% Setting parameters and functionality toggles

% User specify retrodiction methods/ parameters to run...

STFFR_eff_f_t = 1;   % ST-FFR Retrodict with effective frames on/off
retro_eff_fsz = 20;  % Keep same as 'eff_frame_sz'
if STFFR_eff_f_t==1
    load('faces_fsz20','dists_norm'); dists_eff = dists_norm;
end

ind_frm_t = 0;   % Toggle individual frame retrodiction on/off, FFR
STFFR_ind_t = 1;       % ST-FFR toggle
if STFFR_ind_t==1
    load('faces_ind','dists_norm'); dists_ind = dists_norm;
end

sp_grm_t = 1;       % Toggle Speiret's/ Johnstone's method
eff_frame_sz = 20;  % Min 10 for meaningful results of Graeme's algorithm

simulate_t = 0;     % Simulate imaging data. See Sect 2. Else, see Sect 1.

% Data collection specify...
data_sets = 1;       % Number of different data sets to loop over.
eff_frms = 10;       % Number of effective frames for Graeme's method.
                     %      - the number of frames to image over


%% Training parameters
% Set based on same parameters used for training

nb_max=0.0255;  % 0.0255 works well for lab data from Johnstone
res = 255;
round_to = 4;   % Depends on order of magnitude of nb_max in res
fpf = 2;        % as set for training, for provided trained statistics; 
                % fpf=2 as trained in 'faces_ind'.

%% Sect.1 Load experimental imaging data and set reference image

if simulate_t==0

    load('JohnstoneLabData.mat','kept_shots'); % Request from author
    im_sz1 = size(kept_shots,1);
    im_sz2 = size(kept_shots,2);

    ref_image = sum(kept_shots,3);
    
    % Simple hot pixel finding method
    % 2 std dev ~ 80. 2.5 std dev ~ 90. 3 std dev ~ 100. 3.5 std dev ~ 110
    % no hot pixels > 122
    hp_limit = 100;
    [hot_inds_row, hot_inds_col] = find(sum(kept_shots,3)>hp_limit);
    
    ref_image = Hot_pixels(ref_image,hot_inds_row,hot_inds_col);
    ref_image = ref_image./(max(max(ref_image)));

end

%% Sect.2. Loading an input image for measurement simulation

if simulate_t==1
    
    % User specify max n_bar, input image (and image sizes if desired) 
    nb_max_sim = 0.0255; res = 255;
    round_to = 4; % Depends on order of magnitude of nb_max in res
    im_sz1 = 168;
    im_sz2 = 128;
    
    fname='test_target.png';
    input_im=imread(fname); input_im=rgb2gray(input_im); 
    input_im=imresize(input_im,[im_sz1 im_sz2]);
    n_bars_in=double(input_im); lstep_sim = nb_max_sim/res;
    n_bars_in=n_bars_in.*lstep_sim;
    n_bars_in=round(n_bars_in,round_to);

    ref_image = n_bars_in./(max(max(n_bars_in))); % Set reference for ssim
    
    hot_inds_row = [];
    hot_inds_col = [];
end
%% Defining some arrays and functions and setting parameters

cc_ssim = zeros(eff_frms*eff_frame_sz+1,data_sets);% Photon counting
GJ_method_ssim = zeros(eff_frms+1,data_sets);      % Dr Johnstone/ Speirets
retro_ssim = zeros(eff_frms*eff_frame_sz+1,data_sets);       % ST-FFR
retro_ffr_ssim = zeros(eff_frms*eff_frame_sz+1,data_sets);   % FFR
retro_eff_ssim = zeros(eff_frms*eff_frame_sz/retro_eff_fsz+1,data_sets);% ST-FFR-eff

lstep=nb_max/res;
n_bars = 0:lstep:nb_max; n_bars = round(n_bars,round_to);

num_dets = im_sz1*im_sz2;
eta = 0.5; % Lab data recorded with SPAD array with eta = 0.5
RS = ones(im_sz1,im_sz2);
R0 = 1;
pos_t = 0;   % Don't need Bayret calculating pos
eps = 0.005; % As used by Dr Johnstone

len_prior = round(nb_max/lstep + 1);
lstep=nb_max/res;
num_rand_fs = eff_frame_sz*eff_frms;


%% Data collection loops

for data_set=1:1:data_sets
    
    % Choose frames for data set randomly from data if not simulating
    if simulate_t==0
        rand_frame_inds = randi([1 4948],num_rand_fs,1);
        frames = kept_shots(:,:,rand_frame_inds); % Change for no repeats!!(?)
    end
    
    % Initialise Priors and black base images for current data set
    cc_image = zeros(im_sz1,im_sz2); % click counting image
    exp_lam_image = cc_image;        % Graeme/Speirets' method
    prior_ffr = ones(num_dets,len_prior)./len_prior;
    prior_stffr = prior_ffr;
    prior_stffr_eff = prior_ffr;
    
    %% Collect data for each individual frame and effective frame 
    k=1; % Frame kounter
    frame_idx = 1;
    reteff_count = 1;
    for step = 1:eff_frms   % 1 step for 1 effective frame
        
        % Specify current shots (frames) and effective frame
        if simulate_t==0   % Either from data
            shots = frames(:,:,k:k+eff_frame_sz - 1);
            effective_frame = sum(shots,3);
        else   % Or simulate
            [shots,effective_frame] = Measurement(n_bars_in,eff_frame_sz,fpf,im_sz1,im_sz2);
        end

        % Update k for next step
        k = k + eff_frame_sz;

        
        %% Effective frames
        if sp_grm_t == 1
            % Johnstone's method
            [exp_lam,Mj,Mj0,lam_j0p]=Bayret(eta,eps,effective_frame,RS,R0,nb_max,lstep,pos_t);
            
            % Calculate image avg, set 'hot' pixels normalise and calc SSIM
            exp_lam_image = exp_lam_image + exp_lam;
            exp_lam_image_avg = exp_lam_image./step;
            exp_lam_image_avg = Hot_pixels(exp_lam_image_avg,hot_inds_row,hot_inds_col);
            exp_lam_image_avg = exp_lam_image_avg./(max(max(exp_lam_image_avg)));
            GJ_method_ssim(step+1,data_set) = ssim(exp_lam_image_avg,ref_image);
        end
            
        % If performing retro_eff_frm
        if STFFR_eff_f_t == 1

            k_ret = 1;

            for ret_frmd=1:1:eff_frame_sz/retro_eff_fsz

                if retro_eff_fsz~=eff_frame_sz
                    effective_frame = sum(shots(:,:,k_ret:k_ret+retro_eff_fsz - 1),3);
                    [exp_lam,Mj,Mj0,lam_j0p]=Bayret(eta,eps,effective_frame,RS,R0,nb_max,lstep,pos_t);
                    % Update k_ret for next step
                    k_ret = k_ret + retro_eff_fsz;
                end

                [retro_eff]=STFFR(lstep,nb_max,prior_stffr_eff,Mj,Mj0,dists_eff);
                prior_stffr_eff = retro_eff;

                % Calculate image, set 'hot' pixels, normalise and calc SSIM
                retro_eff_image = retro_eff*n_bars.';
                retro_eff_image=reshape(retro_eff_image,im_sz1,im_sz2);
                retro_eff_image = Hot_pixels(retro_eff_image,hot_inds_row,hot_inds_col);
                retro_eff_image = retro_eff_image./(max(max(retro_eff_image)));
                retro_eff_ssim(1+reteff_count,data_set) = ssim(retro_eff_image,ref_image);

                reteff_count = reteff_count + 1;

            end
            
            
        end

        %% Individual frames
        if ind_frm_t == 1
            % For each data frame in effective frame
            for j=1:1:eff_frame_sz

                % Set current frame
                idv_frame = shots(:,:,j);

                % Add clicks from new frame
                cc_image = cc_image + idv_frame;

                % Get Mj and Mjo counts
                [Mj,Nj,Mj0,Nj0,Mbar]=Bayret_mat_Mj(idv_frame,RS,R0);
                % With ST-FFR Retrodictor
                if STFFR_ind_t==1
                    [retro]=STFFR(lstep,nb_max,prior_stffr,Mj,Mj0,dists_ind);
                    prior_stffr = retro;
                else
                    retro = zeros(num_dets,len_prior);
                end
                % With FFR, idv_frm information input
                idv_frm = idv_frame;
                [retro_ffr] = FFR(lstep,nb_max,prior_ffr,idv_frm,fpf);
                prior_ffr = retro_ffr;
                

                % Calculate images from retrodicted distributions
                retro_image = retro*n_bars.';
                retro_ffr_image = retro_ffr*n_bars.';
                retro_image=reshape(retro_image,im_sz1,im_sz2);
                retro_ffr_image = reshape(retro_ffr_image,im_sz1,im_sz2);

                % Set any 'hot' pixels to average of neighbours
                cc_image_norm = Hot_pixels(cc_image,hot_inds_row,hot_inds_col);
                retro_image = Hot_pixels(retro_image,hot_inds_row,hot_inds_col);
                retro_ffr_image = Hot_pixels(retro_ffr_image,hot_inds_row,hot_inds_col);
                
                % Normalising
                cc_image_norm = cc_image_norm./(max(max(cc_image_norm)));
                retro_image = retro_image./(max(max(retro_image)));
                retro_ffr_image = retro_ffr_image./(max(max(retro_ffr_image)));
                
                
                % Calculate ssim
                cc_ssim(frame_idx+1,data_set) = ssim(cc_image_norm,ref_image);
                retro_ssim(frame_idx+1,data_set) = ssim(retro_image,ref_image);
                retro_ffr_ssim(frame_idx+1,data_set) = ssim(retro_ffr_image,ref_image);
                
                frame_idx = frame_idx + 1;
            end
        end

    end
end
%% Calculate data statistics
% Average trajectories
GJ_method_ssim_avgs = mean(GJ_method_ssim,2);
retro_eff_ssim_avgs = mean(retro_eff_ssim,2);
cc_ssim_avgs = mean(cc_ssim,2);
retro_ssim_avgs = mean(retro_ssim,2);
retro_og_ssim_avgs = mean(retro_ffr_ssim,2);

% Standard deviations
GJ_method_ssim_std = std(GJ_method_ssim,0,2);
retro_eff_ssim_std = std(retro_eff_ssim,0,2);
cc_ssim_std = std(cc_ssim,0,2);
retro_ssim_std = std(retro_ssim,0,2);
retro_og_ssim_std = std(retro_ffr_ssim,0,2);


%% Plot/ collect results
if sp_grm_t == 1 && ind_frm_t == 0
    steps_vec = 0:eff_frame_sz:eff_frms*eff_frame_sz;
    figure(1); plot(steps_vec,GJ_method_ssim_avgs,'-x');
    title('Dr Johnstone method Eff frame sz = '+string(eff_frame_sz));
    xlabel('Frames');
    ylabel('SSIM');
end

if ind_frm_t == 1 && sp_grm_t == 1 && STFFR_eff_f_t==0
    steps_vec = 0:eff_frame_sz:eff_frms*eff_frame_sz;
    idv_frames = 0:1:eff_frms*eff_frame_sz;
    method_names = {'Photon Counting';'S/G.J. method';'Trained Retrodictor';'Frm-by-frm Retrodictor'};
    figure(1);
    errorbar(idv_frames,cc_ssim_avgs,cc_ssim_std./sqrt(data_sets)); hold on;%plot(steps_vec,G_method_ssim_avgs,'-x');
    errorbar(steps_vec,GJ_method_ssim_avgs,GJ_method_ssim_std./sqrt(data_sets),'o');%plot(idv_frames,retro_ssim_avgs,'-x');
    errorbar(idv_frames,retro_ssim_avgs,retro_ssim_std./sqrt(data_sets));
    errorbar(idv_frames,retro_og_ssim_avgs,retro_og_ssim_std./sqrt(data_sets));
    legend(method_names,'Location','northwest');
    title('Structural Index of Retrodictive Methods');
    xlabel('Frames'); ylabel('SSIM');
    
end

if ind_frm_t == 1 && sp_grm_t == 1 && STFFR_eff_f_t == 1
    steps_vec = 0:eff_frame_sz:eff_frms*eff_frame_sz;
    idv_frames = 0:1:eff_frms*eff_frame_sz;
    ret_steps_vec = 0:retro_eff_fsz:eff_frms*eff_frame_sz;
    figure(1); 
    method_names = {'Photon Counting';'SIFR';'ST-FFR';'FFR';'ST-FFR-eff'};
    errorbar(idv_frames,cc_ssim_avgs,cc_ssim_std./sqrt(data_sets),'.'); hold on;%plot(steps_vec,G_method_ssim_avgs,'-x');
    errorbar(steps_vec,GJ_method_ssim_avgs,GJ_method_ssim_std./sqrt(data_sets),'x');%plot(idv_frames,retro_ssim_avgs,'-x');
    errorbar(idv_frames,retro_ssim_avgs,retro_ssim_std./sqrt(data_sets),'.');
    errorbar(idv_frames,retro_og_ssim_avgs,retro_og_ssim_std./sqrt(data_sets),'.');
    errorbar(steps_vec,retro_eff_ssim_avgs,retro_eff_ssim_std./sqrt(data_sets),'o');
    legend(method_names,'Location','southeast');
    title('Structural Index of Retrodictive Methods');
    xlabel('Frames');
    ylabel('SSIM');
    
end

if ind_frm_t == 0 && sp_grm_t == 1 && STFFR_eff_f_t == 1
    steps_vec = 0:eff_frame_sz:eff_frms*eff_frame_sz;
    ret_steps_vec = 0:retro_eff_fsz:eff_frms*eff_frame_sz;
    method_names = {'G m';'Retrodictor Eff frm'};
    figure(1);
    plot(steps_vec,GJ_method_ssim_avgs,'-x'); hold on;
    plot(ret_steps_vec,retro_eff_ssim_avgs,'-o')
    legend(method_names,'Location','northwest');
    title('Structural Index of Retrodictive Methods');
    xlabel('Frames');
    ylabel('SSIM');
    
end

if ind_frm_t==1
    if STFFR_ind_t==1
        figure(3); imagesc(retro_image,[0 max(max(retro_image))]); colormap(gray);
        title('Retro Image'); axis off;
    end
    figure(4); imagesc(cc_image_norm,[0 max(max(cc_image_norm))]); colormap(gray);
    title('Click Count Image'); axis off;
    figure(5); imagesc(retro_ffr_image,[0 max(max(retro_ffr_image))]); colormap(gray);
    title('Retro Original Image'); axis off;
end
if sp_grm_t==1
    figure(6); imagesc(exp_lam_image_avg,[0 max(max(exp_lam_image_avg))]); colormap(gray);
    title('Graeme Image'); axis off;
end
if STFFR_eff_f_t==1
    figure(7); imagesc(retro_eff_image,[0 max(max(retro_eff_image))]); colormap(gray);
    title('RetroMLS eff frame Image'); axis off;
end
figure(8); imagesc(ref_image,[0 1]); colormap(gray);
title('Ref Image'); axis off;
