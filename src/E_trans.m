function STERAN = E_trans(Monte, Colli, E)

%     Calculating the collisional step size based on the probability of
%     distribution
%
% 	

Energy = 1:Monte.Step:Monte.Emax;

Step = Monte.Step;
Ediff = 1000;
error = 1.e-5;

%deni = interp1(Energy,Monte.density_states,E/Monte.Step);			% DENI = Density at E
deni = Monte.density_states(round(E/Monte.Step));	
PROB = PDOWN(E,E,Monte.DC, Monte.ICOLL, Monte.Temp); % Evaluate to get initial SCALE
     
cnorme = interp1(Energy,Colli.Cnormal(Monte.ICOLL,:),E/Monte.Step);		%CNORME = Normalization at E
up= interp1(Energy,Colli.Coll_Step(Monte.ICOLL,:),E/Monte.Step);

Tlast = PROB;				% "old" TERM
SUM = 0;
Test = 1;
EE = E;
r_no = rand(1);

if (Monte.IFLAG == 0) 			% Integrate over down steps
      
%******************** Down steps ***************

 ADOWN = r_no*cnorme*(1 - up);	%  Random number-selected integral of down-steps

    while SUM < ADOWN && EE  > 0 && Test > error
             
        H = (EE - Step * round(EE/Step));  	% Align to grain size
        if H <= 0
            H = Step;
        elseif H > EE
            H = EE;
        end
        EE = EE - H;
        %denj = interp1(Energy,Monte.density_states,EE/Monte.Step);
        if EE <= 0
            EE = Monte.Step;
        end
        denj = Monte.density_states(round(EE/Monte.Step));	

        if denj > 0.001
           PROB = PDOWN(E,EE,Monte.DC, Monte.ICOLL, Monte.Temp); % Evaluate to get initial SCALE
           TERM = PROB;
        else
           TERM = 0;
        end

        SPAN = 0.5 * H*(TERM + Tlast);

        if (SUM+SPAN) > ADOWN
            H = 1.01* (ADOWN - SUM)/(0.5* (TERM + Tlast));
            if H < Step
                H = Step;
            end 
            SPAN = 0.5 * H* (TERM + Tlast);  
        end
        SUM = SUM + SPAN;			% Down-step normalization integral
        Tlast = TERM;				% Old value for TERM

        if Tlast > 0.0 && abs(E-EE) > Ediff
           Test = abs(SPAN/SUM);
        end
        
    end

    if Tlast == 0
       EE = EE + H;
    end

        

else						% Up-collisions

% ******************** Up steps ***************

   ETOP = E + 10.*Monte.Temp/1.4388;			%  E + 10 kT
   AUP =  r_no*cnorme*up;			%  Random number-selected

   while SUM < AUP && EE<=ETOP && Test > error

      H = (EE - Step*round(EE/Step));  	% Align to 25 cm-1 grain
      if H <= 0
         H = Step;
      end 
      
      EE = EE + H; 
      %denj = interp1(Energy,Monte.density_states,EE/Monte.Step);
              denj = Monte.density_states(round(EE/Monte.Step));	
    
      if denj >= 0.001
          PROB = PDOWN(EE,E,Monte.DC, Monte.ICOLL, Monte.Temp); % Evaluate to get initial SCALE
          B = exp(-(EE-E)*1.4388/Monte.Temp);			% EE > E
          RATIO = (denj/deni);
          TERM = B*PROB*RATIO;
      else
         TERM = 0;
      end
          
      SPAN = 0.5*H*(TERM + Tlast);

      if (SUM+SPAN) > AUP 
         H = 1.01* (AUP - SUM)/(0.5*(TERM + Tlast));
         if ( H < Step) 
             H = Step;
         end
            SPAN = 0.5*H*(TERM + Tlast);
      end
          
      SUM = SUM + SPAN;			% Up-step normalization integral
      Tlast = TERM;			% Old value for TERM
      if (Tlast > 0 && abs(EE-E) > Ediff) 
          Test = abs(SPAN/SUM);
      end
          
     
  end
       
if (Tlast == 0) 
   EE = EE - H;
end
if EE < 0
   EE = 0;
end

end

STERAN = abs(E-EE);  
      
end
