function [TAU, P] = Stochastic_calc(Monte,Microwave, Collisional, E)

  
%% Calculating reaction rates, timestep and probability of each reaction
%
% FIND TAU 


EC = 40000;
sigma_0 = Microwave.epsi*1e8;
mw_rate = Microwave.MW_power/15.92*sigma_0*1e6;

Energy = 1:Monte.Step:Monte.Emax;

% Calculate reaction rate based on the current temperature


 Rate_1 = 4*Collisional.Rate_1;
 Rate_2 = 4*Collisional.Rate_2;
 RRKM = Collisional.RRKM;
%RKE = interp1(x,v,xq);
        
P(1) = 0.0;   
P(2) = 0.0;     
PUMPU = 0 ;
PUMPD = 0;

if Microwave.MW_on == 1    % Microwave 
     
     if E <= 0
            E = Monte.Step;
     end
     %rho  = Monte.density_states(round(E + Monte.Step)/Monte.Step);
     rho  = interp1(Energy,Monte.density_states,(E +  10*Microwave.MW_Freq));
     rho2  = interp1(Energy,Monte.density_states,E );

     %rho2 = Monte.density_states(round(E)/Monte.Step);
     Ratio = rho2/rho;
     PUMPU = mw_rate; 
     PUMPD = Ratio*PUMPU;
      
end

TARGET = log(1/rand(1));   % for gillespies

if E < EC
    RRKM_E = 0;
else
    RRKM_E =  RRKM(round(E/Monte.Step));
end


SUM = Rate_1 + Rate_2 + RRKM_E;    % Total rate (no pump)
SUM1 = SUM + PUMPU + PUMPD;
TAU = TARGET/SUM1;

%TAU = TARGET/SUM; 

     P(1) = PUMPU*TAU;
     P(2) = PUMPD*TAU;
     P(3) = Rate_1*TAU;        % Collision type 1
     P(4) = RRKM_E*TAU ;       % Unimolecular rxn
     P(5) = Rate_2*TAU;        % Collision type 2
     SUM = sum(P);
     P(1) = P(1)/SUM;

     for i = 2:5
       P(i) = P(i-1) + P(i)/SUM;    % P's compared later with random
     end                           %   numbers: to choose channel


end

