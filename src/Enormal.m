function [Collisional] = Enormal(Monte, Collisional)

% Enormal: Calculating the normalized probability distribution of the collisional
% system

T = Monte.Temp; % Temperature in Kelvin
Egrain = Monte.Step;
iter = 10;
Emax = Monte.Emax;
E = 0:Egrain:Emax-1;
Nmax = Monte.Nmax;
%[Pdist] = Prob_dist(E,Eint,Egrain);
Cup = zeros(Nmax,1);
Cn = ones(1,Nmax);
%Collisional.Coll_Step(Monte.ICOLL,:) = zeros(2,1);
%Collisional.Cnormal(Monte.ICOLL,:) = zeros(Nmax,1);

%E_space= Egrain*i:Egrain:Eint;


% Down steps
for i = 2:Nmax    

    
    E_integrate = 0:Egrain:E(i);

    dens_Ei = Monte.density_states(round(E(i)/Egrain));
    if (dens_Ei >0.01)
        Pdist = PDOWN(E(i),E_integrate,Monte.DC, Monte.ICOLL, Monte.Temp);
        Cup(i) = trapz(E_integrate, Pdist);   
    end
    %Pdist = exp(-(E(i)-E_integrate)./(C1 + C2*E_integrate));
    
    
end

%Shiting all energies by KT  

Cnormal = Cup + T/1.4388;

for j = 1:iter
    
    for i = 2:Nmax
        
        Cn(i) = Cup(i);
        cnorm(i) = Cnormal(i);
        dens_Ei = Monte.density_states(round(E(i)/Egrain));
        
        if (dens_Ei>0.01)
            ETOP = E(i) + 10*T/1.4388;
            EE = E(i);
            if EE < Emax
                dens_E = Monte.density_states(round(EE/Egrain));
                if (dens_E> 0.01)
                    cnormee(i) = Cnormal(round(EE/Egrain));		% Normalization at EE
        
                    ratio = cnorm(i)*dens_E/(dens_Ei*cnormee(i));
                    E_integrate = E(i):Egrain:ETOP;

                    %Pdist = exp(-(E_integrate- EE)./(C1 + C2*E_integrate)); 
                    Pdist = PDOWN(E_integrate, EE,Monte.DC, Monte.ICOLL, Monte.Temp);
                    MulT = Pdist.*exp(-(E_integrate-EE)*1.4388/T);
                    
                    Cn_trap = trapz(E_integrate,MulT);
                    
                    Cn(i) = Cn(i) + Cn_trap;
                else 
                    Cn(i) = Cn(i) + 0;
                
                end
            end
        end

    end
    
    Cnormal = Cn;
end


Coll_step = 1 - Cup./Cn';
   
Coll_step(1) = 0; 

Collisional.Coll_Step(Monte.ICOLL,:) = Coll_step;
Collisional.Cnormal(Monte.ICOLL,:) = Cnormal;

hold on
plot(Coll_step)

 end

