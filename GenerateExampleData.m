%Generate Data for the exponential/gaussian kernel implementation
%Belen Franzoni
%April 14, 2023
%Generates two puntual T1 and T2 populations with noise
clc;clear all; close all;

T11=0.1;%normpdf(tau1,10,1);%
T12=50;

T21=0.5;
T22=20;

d1=100;
d2=10000;      % number of bins in relaxation time grids


tau1 = logspace(-3,4,d1);  
tau2 = logspace(-3,3,d2);  
tau1=tau1';
tau2=tau2';
M=zeros(d1,d2);

for k=1:d1    
    for i=1:d2     
                R11(k)=(1-exp(-tau1(k)/T11));
                R12(k)=1-(exp(-tau1(k)/T12));
                R21(i)=exp(-tau2(i)^2*(1/T21)^2);
                R22(i)=exp(-tau2(i)/T22);
           
    end
end
R11=R11';
M1=R11*R21;
R12=R12';
M2=R12*R22;
M=M1+M2;
figure;
plot(tau2,M(end,:));
figure;
plot(tau1,M(:,1));
figure(6)
surf(tau2,tau1,M)

maxData=max(max(M));
noise=0.05*maxData;
for jj=1:d2
      for ii=1:d1
            Mnoi(ii,jj)=M(ii,jj)+noise*(2*rand-1);
        end
end

save('tauind.dat','tau1','-ascii');
save('tauDir.dat','tau2','-ascii');
save('DataNoise.dat','Mnoi','-ascii');
