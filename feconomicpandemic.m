%%  Chang, Martinez B and Velasco (2020)
% This function solves the economic model 
% It takes as given phi and ht

function [VF]=feconomicpandemic(parameters, SIMUL, wt,et,kt, VFss, sigma)
options = optimoptions('fmincon','Display','off');

%%% Input
% SIMUL is the matrix with elements from the Epidemic Model (phi, ht, p)
% w is the wage vector
% e is the endowment vector
% k is the k vector
% parameters is a matrix of parameters (T kap q D rho cont betta c delta M)
% Matrix VFss with the ss value functions VFss=[vzss vhss vqss vdss vsss];



%%% Output
% Value Functions VF=[vz vh vq vd vs phat p1 p0 Netgain InfRisk MCost MCostInfRisk];
%p1 is the FOC with pt=1
%p0 is the FOC with pt=0
%  Netgain is the current marginal utility net gain of market activities
% InfRisk is the net infection risk
% MCost is the marginal expected Net utility cost of  market activities
% (future)
% McostInfRisk is Mcost*InfRisk (rhs of derivative of maximand)



%% Parameters%%%%%%
    T     =parameters(1,1); % length of simulation
    kap   =parameters(1,2); % kappa
    q     =parameters(1,3); % q
    D     =parameters(1,4); % number of days in hospital
    rhow  =parameters(1,5); % rho work
    contw =parameters(1,6); % gamma work
    betta =parameters(1,7); % Just days, discount factor should be one
    c     =parameters(1,8); % ratio betweeen w and e if both are constant
    d     =parameters(1,9); % Probability of recovering after D days at hospital
    M     =parameters(1,10); % Deadweight loss of death
    rhoh  =parameters(1,11); % rho work
    conth =parameters(1,12); % gamma work
    
    
    
    
%% Parameters (Vectors)
	w=wt; %wage sequence
	e=et; % home endowment
	kappa=kt;  % Kappa vector
	delta=ones(T+D+1,1).*d; % Probability of surviving vector

% Information from  SIAR Model
    phiw=SIMUL(:,8); % Probability of infection - work
	phih=SIMUL(:,9); % Probability of infection - home
	ht=SIMUL(:,1);   % Healthy
	pt=SIMUL(:,10);  % pt
    
%% Empty Vectors for Value Functions
 
vz=zeros(T+D,1);
vh=zeros(T,1);
vq=zeros(T,1);
vd=zeros(T,1);
vs=zeros(T,1);



%% Empty Vectors for other Measures
phat=zeros(T,1);
p1=zeros(T,1);
p0=zeros(T,1);

%% Initial Data (Backward looking)
for j=0:D
vz(T+j,1)=VFss(1,1);
end
vh(T,1)=VFss(1,2);
vq(T,1)=VFss(1,3);
vd(T,1)=VFss(1,4);
vs(T,1)=VFss(1,5);

phat(T,1)=1;


%% Current Consumption  Hospitalized
for i=1:D
    vbetta(i,1)=betta^(i-1);
end

temp=ones(D,1).*e(T,1);
E=[e;    temp];

for i=1:T-1
up=i+D-1;
A=E(i:up,1);
    if sigma==1
        cuh=log(A);
    else
        cuh=A.^(1-sigma)./(1-sigma);
    end
C=cuh.*vbetta;
F(i,1)=sum(C);
end




%% Main recursion
for i=1:T-1
 
   % Current Utility Market
     if sigma==1
     cum=log(w(T-i,1));
     else
     cum=(w(T-i,1).^(1-sigma))./(1-sigma);
     end
     % Probabilities
     probh=kappa(T-i,1).*(1-ht(T-i,1))+kappa(T-i,1).*ht(T-i,1).*phiw(T-i,1); % For q
     
 %Recovered
 vz(T-i,1)=cum+betta*vz(T-i+1,1);    
 
 % Hospitalized
 vh(T-i,1)=F(T-i,1)+betta^(D).*delta(T-i+D,1)*vz(T-i+D,1)-betta^(D).*(1-delta(T-i+D,1))*M;
  
 %Q- Essentials
 vq(T-i,1)=cum+betta.*(probh.*vh(T-i+1,1)+(1-probh).*vs(T-i+1,1));
 
 %Decision Makers
   %Optimal Choice of working
        [p1(T-i,1)]=fBellmanCMVFOC(1,parameters,kappa(T-i,1), ht(T-i,1),w(T-i,1),e(T-i,1),phiw(T-i,1),phih(T-i,1),vh(T-i+1,1),vs(T-i+1,1),sigma);
        [p0(T-i,1)]=fBellmanCMVFOC(0,parameters,kappa(T-i,1), ht(T-i,1),w(T-i,1),e(T-i,1),phiw(T-i,1),phih(T-i,1),vh(T-i+1,1),vs(T-i+1,1),sigma);
        
        if  p1(T-i,1)>0  %Boundary p=1
            phat(T-i,1)=1;  
        elseif p0(T-i,1)<0  %boundary p=0
            phat(T-i,1)=0;
        else  % Interior Solution
            
        fun=@(p)fBellmanCMVFOC(p,parameters,kappa(T-i,1), ht(T-i,1),w(T-i,1),e(T-i,1),phiw(T-i,1),phih(T-i,1),vh(T-i+1,1),vs(T-i+1,1),sigma);
        x0=1;
        phat(T-i,1) =fzero(fun,x0); % Find pt that sets to zero the Bellman's FOC
        
        end
        
       vd(T-i,1)=fBellmanCMV(phat(T-i,1),parameters,kappa(T-i,1), ht(T-i,1),w(T-i,1),e(T-i,1),phiw(T-i,1),phih(T-i,1),vh(T-i+1,1),vs(T-i+1,1),sigma);
  

 % Vulnerable
  vs(T-i,1)=q*vq(T-i,1)+(1-q)*vd(T-i,1);
  
end

%%% Incentives %%%
Netgain=(w-e).*w.^(-sigma);
InfRisk=SIMUL(:,1).*(SIMUL(:,8)-SIMUL(:,9));
vst1=vs(2:T,:);
vht1=vh(2:T,:);
MCost=betta.*kappa(1:T-1).*(vst1-vht1);
MCostInfRisk=MCost.*InfRisk(1:T-1);
p1(T,1)=p1(T-1,1);
p0(T,1)=p0(T-1,1);
MCost(T,1)=MCost(T-1,1);
MCostInfRisk(T,1)=MCostInfRisk(T-1,1);
   

%%%%% Output  %%%%%%%%%%
VF=[vz(1:T,1) vh vq vd vs phat p1 p0 Netgain InfRisk MCost MCostInfRisk];



end