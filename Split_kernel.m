% Fast 2D relaxation distribution estimation with Gaussian and exponential 
%kernel in the transverse relaxation domain.
%
% Adapted from https://github.com/paultnz/flint
% If you use this software, please cite P.D. Teal and C. Eccles. Adaptive
% truncation of matrix decompositions and efficient estimation of NMR
% relaxation distributions. Inverse Problems, 31(4):045010, April
% 2015. http://dx.doi.org/10.1088/0266-5611/31/4/045010
% 
% AND the implementation presented by M.B. Franzoni, G.S. Vila , E.V. Sillettaa,  G.A. Monti , D. Masiero , R.H. Acosta, E.A. Domenéc 

t_c=0.6;   % cuttoff time between Gaussian and exponential kernels

alpha=5E0;  % Regularization parameter

Nx = 100;   % number of bins in relaxation time grids
Ny = 100;

T1 = logspace(-3,4,Ny);  %T1  inverted domain
T2 = logspace(-3,3,Nx);  %T2  inverted domain

data=load('DataNoise.dat');
tau1=load('tauind.dat');
tau2=load('tauDir.dat');


% Saturation-Recovery kernels
K1=  1-exp(-tau1*(1./T1));

Z=data;
% Transverse relaxation kernels
for tt=1:size(tau2)
    for uu=1:Nx
        if T2(uu)<t_c
            K2(tt,uu) = exp(-tau2(tt)^2 *(1/T2(uu))^2 );
        else
            K2(tt,uu) = exp(-tau2(tt)*(1/T2(uu))); 
        end
    end
end


% 2D inversion may be downloaded from:  https://github.com/paultnz/flint

[S,resida] = flint(K1,K2,Z,alpha);%S is the inverted relaxation matrix.

figure(10)
contour(T2,T1,S,90)
set(gca,'YScale','log','FontSize',13)
set(gca,'XScale','log','FontSize',13)
xlabel('T_{2} [ms]','FontSize',18)
ylabel('T_{1} [ms]','FontSize',18)
caxis([0 0.3])
colorbar


