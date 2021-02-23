%% Chang, Martinez and Velasco (2020)
% Equation 8 - Bellman Equation

function [f]=fBellmanCMVFOC(p,parameters,k, h,w,e,phiw,phih,vht1,vst1,sigma)

% p
% parameters (betta)
%  k - probability of getting symptoms
% h   - share of healthy
% w   - market reward
% e   - home endowment
% phiw - infection risk at market
% phih- infection risk at home
% vh  - Value function of hospital one period ahead
% vs  - Valiue fuction of vulnerable one period ahead
% sigma - Preferences of individuals

%parameters
betta=parameters(1,7);
kappa=k;

% consumption
c=p.*w+(1-p).*e;

% Current Utility
    if sigma==1
        cufoc=1./c;
    else
        cufoc=c.^(-sigma);
    end



% FOC Bellman Equation
f=cufoc.*(w-e) - betta.*kappa.*h.*(phiw-phih).*(vst1-vht1);










end
