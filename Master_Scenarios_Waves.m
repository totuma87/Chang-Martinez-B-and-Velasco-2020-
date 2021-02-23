% Chang, Martinez and Velasco (2020)
% Base Scenario

clc
clear
 
% Base Calibration Parameters 
run CalibrationBase.m;

% Changes in  Calibration Parameters
iteration=1;
for h=[150]
    
    cmax=0.99; %https://www.cnbc.com/2020/05/09/it-pays-to-stay-unemployed-that-might-be-a-good-thing.html
         %https://fivethirtyeight.com/features/many-americans-are-getting-more-money-from-unemployment-than-they-were-from-their-jobs/#:~:text=The%20idea%20behind%20a%20%24600,recipients%20was%20%24970%20per%20week.
	cmin=0.10;
    
    
     for j=30:h
       cbar(j,1)=cmax;
     end
     e=w.*cbar; % home endowment
    
  
  %% Solving Decentrilized Model
 % Covergence and Model Criteria
   maxiter=700000; % Maximum number of Iterations to try to reach equilibrium
   lambda=0.99995; % measure for next guess % This item turned out to be very important. We want a smooth change between pt iterations
   concriteria=0.01; % Convergence Criteria
   
   % Initial Guess
   if iteration==1
   load('ScenarioBase.mat','SIRF')
   p0=SIRF(:,10); 
   else
   load('p0INT.mat','poINT')
   p0=poINT;
   end
   
   % Shooting Algorithm
  [SIRF,VFssF, VFF, flag, ~, ~]=fdecentralizedshooting(p0, w, e, kappa,parameters,sigma, initialvalues,lambda,concriteria,maxiter)
  poINT=VFF(:,6);
  save('P0INT.m','poINT')

toc
  
  
  
  
iteration=iteration+1;  
  
  
  
  %save(['Scenario_Waves' num2str(h) '.mat'],'SIRNoE','SIRNOP','SIRDC') 
  save(['Scenario_Waves_m' num2str(h) '.mat'],'SIRF', 'VFF') 

end


