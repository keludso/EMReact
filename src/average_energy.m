function [Avg_energy] = average_energy(edist, Energy)

%average_energy: Plots the average energy in the simulation interval 
%   edist: Energy distribution matrix - [NMAX,TLIM]
%   Energy: 0:Egrain:EMAX


Avg_en = 0
for tm = 1:101
    Avg_en = 0
    for enm = 1:8000

        Avg_en = Avg_en + edist(tm,enm)*Energy(1,enm);
        
    end

    Avg_energy(tm) = Avg_en/2000;
end

%figure
%plot(Avg_energy)
end

