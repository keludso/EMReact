function [status,E, NET_A, NET_B, Monte] = monte_calc(Monte, Colli, E, P,tau, Microwave)


%% monte_calc: This function calculates the reactiuon that will occur based on the
% probability of the reaction using he Gillespies algorithm 

if (E< Monte.Step)
    E = Monte.Step;
end


Monte.ICOLL = 1;
NET_A = 0;
NET_B = 0;
MW_absorbed =0;
RX = rand(1);
Energy = 1:Monte.Step:Monte.Emax;

% Glycerol
loss_tangent = 0.540;
Cavity_area = 0.01;
Volume = 22.4;
density = 1.1;
ang_omega = 0.7;

%field = Microwave.MW_power/(4*pi*(3e10/Microwave.MW_Freq)*Microwave.epsi*Volume);
field = sqrt(2*Microwave.MW_power/(Cavity_area*Microwave.epsi0*3e8));
dielectric_relax = 1e-8;

Power_dissipated = (2.4e9*Microwave.epsi*loss_tangent)*(field^2)*5.2e22*tau/6.23e23;

microwave_heat = (2*pi*3e10/Microwave.MW_Freq)*Microwave.epsi*Microwave.loss_tangent*field/(density*Colli.CVI);

if Microwave.MW_on > 0
       if RX < P(1)                   % Photon absorption
            MW_jump = Power_dissipated;  %cm^-1
            E = E + MW_jump; %cm^-1
            NET_A = 1;
            status = 0; % Back to the start  
            return;
    
       elseif RX < P(2)        % Stimulated emission
           MW_jump = 0.3*Power_dissipated;
           E = E - MW_jump ;               % Photon stimulated emission
           if E < 0 
              E = 0;
           end
           NET_B = -1;
           status =0;  % Back to the start 
           return;
       end
            
end

      if RX < P(3)           % Collision type 1
          %RNTER = interp1(Energy,Colli.Coll_Step(Monte.ICOLL,:),E/Monte.Step);
          RNTER = Colli.Coll_Step(Monte.ICOLL,round(E/Monte.Step));
          RXY = rand(1);
          if (RXY < RNTER || E <= 0)    	% An up-step was chosen
              Monte.IFLAG = 1;
             STERAR = E_trans(Monte, Colli, E);
             E = E + STERAR;
          else
            Monte.IFLAG = 0;  
            STERAR = E_trans(Monte, Colli, E);
            E = E - STERAR;
          end
          Monte.Temp = Monte.Temp + microwave_heat*tau;
         status = 0;
         return;
          
      elseif RX < P(4)           % Unimolecular rxn
          
            
            status =1;
            return;
      
      elseif RX < P(5)           % Collision type 2
          Monte.ICOLL = 2;
          %RNTER = interp1(Energy,Colli.Coll_Step(Monte.ICOLL,:),E/Monte.Step);
          RNTER = Colli.Coll_Step(Monte.ICOLL,round(E/Monte.Step));
          RXY = rand(1);
          if (RXY < RNTER || E <= 0)    	% An up-step was chosen
             Monte.IFLAG = 1;  
             STERAR = E_trans(Monte, Colli, E);
             E = E + STERAR;
          else
            Monte.IFLAG = 0;    
            STERAR = E_trans(Monte, Colli, E);
            E = E - STERAR;
          end
         Monte.Temp = Monte.Temp + microwave_heat*tau;
         status = 0;
         return;           
      end
     
    
status =0;

end


