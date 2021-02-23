%% Chang, Martinez and Velasco (2020)
% Equation 8 - Bellman Equation

function [f]=fBellmanCMV(p,parameters,k, h,w,e,phiw,phih,vht1,vst1,sigma)

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
        cu=log(c);
    else
        cu=(c.^(1-sigma))./(1-sigma);
    end
% Phibar
phi=p.*phiw +(1-p).*phih;

% probability of going to hospital tomorrow
probh=kappa.*(1-h)+kappa.*h.*phi;

% Bellman Equation
f=cu + betta.*(probh.*vht1+(1-probh).*vst1);










end
