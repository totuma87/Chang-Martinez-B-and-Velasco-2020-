%% Chang, Martinez, and Velasco (2020)

function [f]=feconomicOP(pt,parameters,w,e,initialvalues, sigma)

% Parameters

T=parameters(1,1);       % Length of Simulation
k=parameters(1,2);       % Probability of getting symptoms
q=parameters(1,3);       % Share of Essentials
D=parameters(1,4);       % Days Hospitalized
rhom=parameters(1,5);    % Number of contacts  market
gammam=parameters(1,6);  % Transmission rate market
betta=parameters(1,7);   % Betta
delta=parameters(1,9);   % Probability of Survival
M=parameters(1,10);   % Cost of Life
rhoh=parameters(1,11);   % Number of Contacts household
gammah=parameters(1,12); % Transmision rate household
Ms=parameters(1,13);   % Social Cost of life


% SIR Model 
[SIR]=fpandemic(pt,parameters,initialvalues);

% SIR Model 
h=SIR(:,1);
s=SIR(:,2);
H=SIR(:,3);
x=SIR(:,4);
z=SIR(:,5);
J=SIR(:,6);
omega=SIR(:,7);

% Steady State
	% Assumption pss=1
	if sigma==1
	css=log(w(T,1));
	else
	css=(w(T,1).^(1-sigma))./(1-sigma);
	end

	VFss=(1-omega(T,1)).*css./(1-betta);

% Vector Empyt
VF=zeros(T,1);

% Initial Values
VF(T,1)=VFss;


% Main Recursion
for t=1:T-1

	c=pt(T-t,1)*w(T-t,1)+(1-pt(T-t,1))*e(T-t,1);

		if sigma==1
				cu=log(c);
				wu=log(w(T-t,1));
				eu=log(e(T-t,1));
		else
				cu=(c.^(1-sigma))./(1-sigma);
				wu=(w(T-t,1).^(1-sigma))./(1-sigma);
				eu=(e(T-t,1).^(1-sigma))./(1-sigma);
	    end

	U=s(T-t,1)*q.*wu+z(T-t,1).*wu-J(T-t,1).*Ms+x(T-t,1).*eu+s(T-t,1).*(1-q).*cu;

	VF(T-t,1)=U+betta.*VF(T-t+1,1);

end
f=VF(1,1);
end







  



