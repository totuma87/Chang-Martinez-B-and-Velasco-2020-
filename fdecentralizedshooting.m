% Chag, Martinez B, and Velasco (2020)
% Shooting Algorithm Following Pissarides et al (2020)

% Input
% p0 Initial guess for pt
% w  Market reward vector
% e Home reward vector e 
% kappa Vector of kappa
% parameters - Parameters parameters=[T kap q D rhow contw betta c delta M rhoh conth Ms];
% sigma  _ Preferences
% initialvalues - Initial epidemiological Values initialvalues=[h0 s0 H0 x0 z0 J0 omega0 xi0];
% lambdap  - Smoothing parameter for pt between iterations
% concriteriap - Convergence criteria for euclidean distance between pt
% maxiterp - Max amount of iterations to find equilibrium

% Output
% SIRF Matrix with epidemiological results and 
% VFssF Row vector with steady state value functios
% VFF   Matrix with Value Functions
% flag   % Result of algorithm 1 converged, 0 ran out of iterations
% piter   % pt of iterations
% convergence % Euclidean distance of iterations

function [SIRF,VFssF, VFF, flag, piter, convergence]=fdecentralizedshooting(p0, w, e, kappa,parameters,sigma, initialvalues,lambdap,concriteriap,maxiterp)
% Covergence and Model Criteria
   maxiter=maxiterp; % Maximum number of Iterations to try to reach equilibrium
   lambda=lambdap; % measure for next guess % This item turned out to be very important. We want a smooth change between pt iterations
   concriteria=concriteriap; % Convergence Criteria
   T=parameters(1,1);
   
  % Empty Matrix to hold econonmic pt of iterations
    %piter=zeros(T,maxiter);
    %ptrans=zeros(T,maxiter);
    %convergence=zeros(1,maxiter);
   
   gamma=0;
   iter=1;
   iter2=1;
   loop=0;
   
   % Initial Guess for p(t)
   p=p0;  % p(t)
       

 while gamma<1 && iter2<=2
           loop=1+loop;
           
          
           
           'Iteration'
           loop
                % Epidemic Model
               %[SIR, SIMUL]=fsimpandemic(parameters,kappa, initialvalues, p, dummyrecov);
               [SIR]=fpandemic(p,parameters,initialvalues);

                % Steady States
               [VFss]=fdynamiceconomicss(SIR, w, e, kappa, parameters, sigma);

                % Model
               [VF]=feconomicpandemic(parameters, SIR, w,e,kappa, VFss, sigma);
               
               ptrans(:,loop)=p;
               piter(:,loop)=VF(:,6);
               distance=sqrt(sum((p-VF(:,6)).^2,1));
               convergence(:,loop)=distance;
                               
               'Convergence'
                    distance
                        
       if distance<=concriteria

                   'The Model reached an equilibrium'
                     gamma=1;
                     flag=1;       
                     [SIRF]=fpandemic(VF(:,6),parameters,initialvalues);
                     [VFssF]=fdynamiceconomicss(SIRF, w, e, kappa, parameters, sigma);
                     [VFF]=feconomicpandemic(parameters, SIRF, w,e,kappa, VFssF, sigma);

                                          
                     


               else
                  'One more Iteration'
                   p(:,1)=lambda*p(:,1)+(1-lambda)*VF(:,6);
                              
                  
                  
                  iter=iter+1;
                  if iter>maxiter
                     'Run out of Iteration attempts without reaching equilibrium - Lets try new initial guess for pt'
                   iter2=3;
                   iter=1;
                      if iter2>2
                           'Run out of Iteration attempts without reaching equilibrium'
                           flag=0;
                      end 
                   end
     end
 end
 
 end
