%%  Chang and Velasco (2020)
% By Humberto Martinez B - Economic Model
% This function calculates steady states of value functions
% It is assumed that in the ss, Vm>Ve (pss=1)
%  SS values of wages, endowments, kappa, and SIR variables are in the last
%  rows of vectors
% This version includes risk aversion function

function [VFss] =fdynamiceconomicss(SIR, w, e, k, parameters, sigma)



%%% Input
% SIR is the matrix with elements from the Epidemic Model (phi, ht, pw,ph)
% w is the wage vector
% e is the endowment vector
% k is the k vector
% parameters is a matrix of parameters (T kap q D rhow contw betta c delta M rhoh conth)

%%% Output
% Matrix VFss with the ss value functions [vrss vhss vmss vess
% vdss vsss]


% Parameters
T=parameters(1,1);
D=parameters(1,4);
betta=parameters(1,7);
delta=parameters(1,9);
q=parameters(1,3);
M=parameters(1,10);

% Picking ss values
%hss=SIR(T,2);
%phiwss=SIR(T,1);
%phihss=SIR(T,4);
hss=SIR(T,1);
phiwss=SIR(T,8);
phihss=SIR(T,9);

wss=w(T,1);
ess=e(T,1);
kss=k(T,1);

%Value Functions
% Assume that p^{ss}=1, phiw==0, hss=1 - This is true for plausible
% parameter values of rho and gamma

% Current utility
if sigma==1
    
    cum=log(wss);
    cuh=log(ess);

else
    cum=(wss.^(1-sigma))./(1-sigma);
    cuh=(ess.^(1-sigma))./(1-sigma);
    
end

% Probabilities
probh=kss.*((1-hss)+hss.*phiwss);

% Value Functions Steady States
vzss=cum./(1-betta);
vsss=cum./(1-betta);
vhss=(cuh.*(1-betta^D))./(1-betta) + betta^D.*(delta.*vzss-(1-delta).*M);
vqss=cum./(1-betta);
vdss=cum./(1-betta);




VFss=[vzss vhss vqss vdss vsss];

end
