function [RKCI1,RKCI2, ALNRKE] = calculate_rate(Monte,Collisional)

% calculate_rate: Calcuating the reaction rates of the system based on the
% reaction parameters provided in the file. 
% Monte , Collisional : Object class defined with inital parameters.


% Reaction rate 1 
RKCI1 = (8.09e-10)*Collisional.PAPRES*(sqrt(Monte.Temp/1000))*(sqrt(40/Collisional.AMASS1))*((Collisional.SIG1/5)^2)/(0.636 + (0.246*log(Monte.Temp/Collisional.WELL1)));
% Reaction RAte 2 
RKCI2 = (8.09e-10)*Collisional.BGP*(sqrt(Monte.Temp/1000))*(sqrt(20/Collisional.RMASS))*((Collisional.XSECT/5)^2)/(0.636+ 0.246*log(Monte.Temp/Collisional.WELL));

EC = 45000;
%% RRKM parameters for UNIMOLECULAR DISSOCIATION

%     IRATE = 0 : INVERSE LAPLACE RRKM, 
%           = 1 : RRKM.  Input 300 element double-array of ln(sum of states) in ALNRKE
%     AVENU = LOG10(AFACTOR)
%     EC = critical energy (cm-1)
%     ROTDGN =  rotation factor*path degeneracy
%
% K(E) = A*density(E-Ecrit)/density(E)

Energy = 1:Monte.Step:Monte.Emax;
ECRIT = Monte.Emax - EC;
AVENU = 17.3;

for j=2:length(Monte.density_states)

       E = (j-1)*Monte.Step;
       %E  = EE(j-1);
       
     

       %DENS = interp1(Energy, Monte.density_states, E/Monte.Step); % Calculate density of states by interpolating the between values 
       DENS = Monte.density_states(round(E/Monte.Step));
       DENOM = log(DENS);
       %DENS = interp1(Energy, Monte.density_states,EE/Monte.Step);    % for Laplace Transform k(E)
       DENS = Monte.density_states(E/Monte.Step);
       ANUM = log(DENS);
       ALNRKE(j) =  2.302585*AVENU + ANUM - DENOM;	% Inverse Laplace method
       %ALNRKE(j) = log(AMASS2*5.2)- 24534/Monte.Temp - (2.3*log(Monte.Temp/1000))  + DENOM

       
end
        

end

