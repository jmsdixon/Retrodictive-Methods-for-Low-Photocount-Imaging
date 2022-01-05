% Author: Dr Graeme Johnstone, University of Strathclyde.

function [exp_lam,Mj,Mj0,lam_j0p]=BayretG_RSmat(eta,eps,pic,RS,R0,lmax,lstep,pos_t)
% This code is based on the work by Speirits et al. It is used to
% reconstruct low photon count iamges of objects illuminated by structured
% light sources. In the case that RS is set to a matrix of ones, it will
% recreate the work of Speirits et al. exactly. It calculates two outputs
% for a given input picture.
% exp_lam, the expectation value of the light intensity at a given pixel
% pos, the posterior probability distribution for the intensity of light at
% that pixel

% the inputs are as follows
% eta = quantum efficiency of the detector, the probability that a photon
% incident upon the detector will result in a detected photon.
% eps = dark count rate for the detector

% pic = the image that is to be reconstructed
% RS = the illumination pattern used to with this image

% R0 = Radius defining the size of the measurement window and therefore the
% size of the window for the prior. The higher this number is, the "smoother"
% the final image will be

% lmax = the maximum intensity that the posterior should be calculated
% for
% lstep = the step size for calculating posterior probability

% pos_t is toggled 1, to calculate the posterior probability or 0 to switch
% it off. This consumes a lot of computational time and I rarely switch it
% on


%% Calculate Mj, Nj, Mj0, Nj0 and mbar

[Mj,Nj,Mj0,Nj0,mbar]=Bayret_mat_Mj(pic,RS,R0);

maxN=max(max(Nj));
Nj=ones(size(Mj)).*maxN;
maxN0=max(max(Nj0));
Nj0=ones(size(Mj)).*maxN0;

% Trying a different mbar, 01/12/20
% mbar=(sum(pic,'All')-(Mj+Mj0))/numel(pic);
% I'm still not happy with mbar, it should be a little less than the
% average value of the image, so I will use that as an approximation
mbar=sum(pic,'All')/numel(pic);

Mj=RS.*Mj;
Mj0=RS.*Mj0;

% calculate the prior from the whole rest of the image lam_j0p
lam_j0p=(mbar-eps)/eta;

%% Calculate some useful intermediate functions

c=Mj0+1; % alpha in the paper
d=Nj0+(1./(eta*lam_j0p)); % beta in the paper

%% Calculate the prior lam_jp, which is dependent on lam_j0p, through c and d

num=gammainc(c+1,eps*d,'upper').*gamma(eps*d);
den=d.*(gammainc(c,eps*d,'upper').*gamma(eps*d));
nod=num./den;
lam_jp=(1/eta)*(nod-eps); % the local prior

%% Calculate some more useful intermediate functions

a=Mj+1;
b=Nj+(1./(eta*lam_jp));

%% Calculate the expectation value for the intensity at pixel j
eb=eps*b;
fix=0.00001;
eb(isnan(eb))=fix; % Definitely cheating, this removes NaNs that are caused by numbers in the denominator being too small
eb(eb<=fix)=fix;
% The previous line is included to correct for an issue caused in gamma
% function of for the denominator (currently at line 43). The problem is
% that when c=~1000, this tends to zero. A zero in the denominator is a
% problem leading to NaNs. There is probably a better correction for this
% but i think this will help for now.
num=((eb).^a).*exp(-a.*b);
den=b.*gammainc(a,eb,'upper').*gamma(eps*b);

exp_lam=(1/eta)*((a./b)+(num./den)-eps);

% Calculate the probability distribution for the intensity at pixel j

pos=zeros(size(pic,1),size(pic,2),round((lmax/lstep)+1));
n=0;
% if eps*b > 0
%     eb=eps*b;
% else
%     eb=0.01;
% end

gam=gammainc(a,eb,'upper').*gamma(eps*b); % This is the denominator of the posterior probabitilty distribution equation and is independent of l, therefore it is calculated initially instead of inside the for loop.

if pos_t==1
    for l=0:lstep:lmax
        n=n+1;
        pos(:,:,n)=Bayret_pos(l,eta,eps,a,b,gam);
    end
end

% This section will be used to correct the data in the dark stripe to zero
% and to set the probabitility distribution to a kronecker delta function
% at zero.

for x=1:size(pic,1)
    for y=1:size(pic,2)
        if RS(x,y)==0
            exp_lam(x,y)=0;
            pos(x,y,:)=zeros((lmax/lstep)+1,1);
            pos(x,y,1)=1;
        end
    end
end