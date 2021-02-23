% Chang, Martinez and Velasco (2020)
% Base Scenario

clc
clear
 
% Base Calibration Parameters 
run CalibrationBase.m;


% Changes in  Calibration Parameters
    % No changes in the base scenario

 %% Epidemic Model No Economy - Save Results
  onesvector=ones(T,1);
  [SIRNoE]=fpandemic(onesvector,parameters,initialvalues);

  %% Solving Decentrilized Model
 
 % Covergence and Model Criteria
   maxiter=300000; % Maximum number of Iterations to try to reach equilibrium
   lambda=0.9995; % measure for next guess % This item turned out to be very important. We want a smooth change between pt iterations
   concriteria=0.01; % Convergence Criteria
   
   % Initial Guess
   p0=onesvector; 
   
   % Shooting Algorithm
  [SIRF,VFssF, VFF, flag, ~, ~]=fdecentralizedshooting(p0, w, e, kappa,parameters,sigma, initialvalues,lambda,concriteria,maxiter)




  
  
  
  
  
  
 %% Optimal Maximization Problem
 options=optimoptions('fmincon','Algorithm','sqp','MaxFunctionEvaluations',300000, 'Display','iter');
 fun=@(p)feconomicOP(p,parameters,w,e,initialvalues, sigma)*(-1);
 p0=SIRF(:,10);
  
 pOP=fmincon(fun,p0,[],[],[],[],zeros(T,1),ones(T,1),[],options);
 [SIROP]=fpandemic(pOP,parameters,initialvalues);
  

save('ScenarioBase','SIRNoE','SIROP', 'SIRF', 'VFF') 





