%% Chang, Martinez, and Velasco (2020)

function [SIR]=fpandemic(pt,parameters,initialvalues)

% Parameters

T=parameters(1,1);    % Length of Simulation
k=parameters(1,2);    % Probability of getting symptoms
q=parameters(1,3);    % Share of Essentials
D=parameters(1,4);    % Days Hospitalized
rhom=parameters(1,5); % Number of contacts  market
gammam=parameters(1,6); % Transmission rate market
delta=parameters(1,9);  % Probability of Survival
rhoh=parameters(1,11);  % Number of Contacts household
gammah=parameters(1,12); % Transmision rate household

% Initial Values 
h0=initialvalues(1,1);
s0=initialvalues(1,2);
H0=initialvalues(1,3);
x0=initialvalues(1,4);
z0=initialvalues(1,5);
J0=initialvalues(1,6);
omega0=initialvalues(1,7);
xi0=initialvalues(1,8:8+D-1);


% Matrices
h=zeros(T,1);
s=zeros(T,1);
z=zeros(T,1);
omega=zeros(T,1);
J=zeros(T,1);
H=zeros(T,1);
phim=zeros(T,1);
phih=zeros(T,1);
x=zeros(T,1);
xi=zeros(T,D);


% Initial Values
h(1,1)=h0;
s(1,1)=s0;
z(1,1)=z0;
omega(1,1)=omega0;
J(1,1)=J0;
H(1,1)=H0;
x(1,1)=x0;
xi(1,1:D)=xi0;


% Main Recursion

for t=1:T-1

%Probabilities
hw=((q+(1-q)*pt(t,1))*H(t,1)+z(t,1))/((q+(1-q)*pt(t,1))*s(t,1)+z(t,1));
phim(t,1)=1-(hw+(1-gammam)*(1-hw))^rhom;
phih(t,1)=1-(h(t,1)+(1-gammah)*(1-h(t,1)))^rhoh;
probbar=(q+(1-q)*pt(t,1))*phim(t,1)+(1-q)*(1-pt(t,1))*phih(t,1);

% Healthy Vulnerables
H(t,1)=s(t,1)*h(t,1);
H(t+1,1)=H(t,1)*(1-probbar);

% Transition Hospitalized
 for j=2:D
	xi(t+1,j)=xi(t,j-1);
 end

% Hospitalized
xi(t+1,1)=k*(s(t,1)-H(t,1)+probbar*H(t,1));
x(t+1,1)=sum(xi(t+1,:),2);

% Vulnerables
s(t+1,1)=s(t,1)- xi(t+1,1);

% Healthy
h(t+1,1)=H(t+1,1)/s(t+1,1);

% Recovered
z(t+1,1)=z(t,1)+delta*xi(t,D);

% Dead
J(t+1,1)=(1-delta)*xi(t,D);
omega(t+1,1)=omega(t,1)+J(t+1,1);
end


% End of Period Probabilities
hw=((q+(1-q)*pt(T,1))*H(T,1)+z(T,1))/((q+(1-q)*pt(T,1))*s(T,1)+z(T,1));
phim(T,1)=1-(hw+(1-gammam)*(1-hw))^rhom;
phih(T,1)=1-(h(T,1)+(1-gammah)*(1-h(T,1)))^rhoh;


% Result Matrix
SIR=[h s H x z J omega phim phih pt xi];

end
  



