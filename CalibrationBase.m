%% Chang, Martinez and Velasco (2020)
%  m-file with Base Calibration


% Parameters
    T = 2500; % length of simulation
    
    sigma=1;  %Risk aversion 1 u(c)=ln(c), >1 u(c)=c^(1-sigma)/(1-sigma) 
    kap = 0.18; % kappa
    kappa=ones(T,1)*kap;  % Kappa vector
    q = 0.3;  % q  - share of essentials
    D = 18;   % number of days in hospital
    
    rhow=9; % WORK - POLYMOD Study (Households 2.52 for USA - Polymod for hh2 10.65 and for HH3 12.87, then for US 11.7 of which 23% were at home)
    contw = 0.05; % WORK - contagion parameter For a R0 of 2.4
    rhoh = 2.7; % HOME - Size of households (Average USA)
    conth = 0.05; % HOME - contagion parameter For a R0 of 2.4
    
    c=0.38; % ratio betweeen e and w if both are constant (https://fivethirtyeight.com/features/many-americans-are-getting-more-money-from-unemployment-than-they-were-from-their-jobs/#:~:text=The%20idea%20behind%20a%20%24600,recipients%20was%20%24970%20per%20week.)
    delta=1-0.01;% Probability of recovering after D days at hospital
    
    
        
    sdf=0.01; % 1% real interest rate
    betta=(1/(1+sdf))^(1/365); 
    % 60062.2 Annual GDP per capita (current $) US WDI 2017
    % $10 million Kneisner and Viscusi (2019)  
    wobj=10000000*(1-betta^365)/60062.2;
    m=1;
    if sigma==1
    cw=log(wobj); % Statistical Value of Life consistent with the model
    else
    cw=(wobj^(1-sigma))./(1-sigma);
    end
    m=1;
    M=cw/(1-betta)*m;
    
    % Social M 
    % Share of Labour Compensation in GDP at Current National Prices for
    % United States (60%)
    ws=wobj/0.6;
    if sigma==1
    cws=log(ws); % Statistical Value of Life consistent with the model
    else
    cws=(ws^(1-sigma))./(1-sigma);
    end
    m=1;
    Ms=cws/(1-betta)*m;
    Ms=M; % For simplicity
    
    
    % Vector of Cbar (the greater c, the better is to stay at home)
    % Average UN benefit is 370,  
    cbar=ones(T,1).*c;
    
    % Rewards Location
    w=ones(T,1); %wage sequence
    e=w.*cbar; % home endowment
    
    
    
    % Initial Values
    h0 = 1-10000/330000000; % Healthy
    s0 = 1;      % Susceptibles
    H0 = s0*h0;
    z0 = 0;      % Recovered
    J0 = 0;        % Deaths at period 0
    omega0 = 0;    % Accumulated Dead
    x0=0;            % Total Hospitalized
    xi0=zeros(1,D); % Hospitalized at period 1 by days hospitalized

    %% Calibration
    parameters=[T kap q D rhow contw betta c delta M rhoh conth Ms];
    initialvalues=[h0 s0 H0 x0 z0 J0 omega0 xi0];
