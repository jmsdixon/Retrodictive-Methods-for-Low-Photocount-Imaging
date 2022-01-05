% Author: Dr Graeme Johnstone, University of Strathclyde.

function [Mj,Nj,Mj0,Nj0,Mbar]=Bayret_matRS_Mj(pic,RS,R0)

R1=R0+1;

M=zeros(size(pic,1)+2*R1,size(pic,2)+2*R1); % M should be the size of the pic + 2*R0
M0=zeros(size(pic,1)+2*R1,size(pic,2)+2*R1); % M0 should be two bigger in both dimensions
N=zeros(size(pic,1)+2*R1,size(pic,2)+2*R1); % N has to be the same size as M
Ni=zeros(size(pic,1)+2*R1,size(pic,2)+2*R1); % I need another matrix to help with the counting!
N0=zeros(size(pic,1)+2*R1,size(pic,2)+2*R1); % N0 has to be the same size as M0

% The pic should be recreated in an intermediate matrix, pici, with a border
% of zeros R1 thick around the outside of the pic
pici=zeros(size(pic,1)+2*R1,size(pic,2)+2*R1);
pici(R1+1:size(pici,1)-R1,R1+1:size(pici,2)-R1)=pic;
% Creating an equivalent matrix for N, here the "real" pixels have a value
% of 1 and the fake pixels have a zero and the pixels in an off stripe also
% have zero.
Ni(R1+1:size(pici,1)-R1,R1+1:size(pici,2)-R1)=1.*RS;

for x=ceil(-R1):floor(R1)
    for y=ceil(-R1):floor(R1)
        r=sqrt(x^2+y^2);
        if r<=R0
            Mint=circshift(pici,[x,y]);
            M=M+Mint;
            Nint=circshift(Ni,[x,y]);
            N=N+Nint;
        elseif r>R0 && r<=R1
            Mint=circshift(pici,[x,y]);
            M0=M0+Mint;
            Nint=circshift(Ni,[x,y]);
            N0=N0+Nint;            
        end
    end
end

Mj=M(R1+1:size(pici,1)-R1,R1+1:size(pici,2)-R1);
Mj0=M0(R1+1:size(pici,1)-R1,R1+1:size(pici,2)-R1);

Nj=N(R1+1:size(pici,1)-R1,R1+1:size(pici,2)-R1);
Nj0=N0(R1+1:size(pici,1)-R1,R1+1:size(pici,2)-R1);

%% A background value for each pixel needs to be calculated

% This value should be equal to the total number of counts in the
% image-(Mj+Mj0)/the total number of pixels-(Nj+Nj0).

Mbar=zeros(size(Mj,1),size(Mj,2));
Mtot=sum(sum(pic));
Ntot=numel(pic);
num=Mtot-(Mj+Mj0);den=Ntot-(Nj+Nj0);
Mbar=num./den;

