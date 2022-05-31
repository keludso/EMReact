classdef Micro_wave
    %Microwave: Object class defining the initial parameters of the
    %microwave parameters.
    
    properties
        MW_on = 0;
        MW_Freq =0;
        MW_power = 0;
        loss_tangent = 0.825;
        epsi0 = 8.854e-12;
        epsi = 8.854e-12*46.7;
        Volume = 0.2^3;
        density = 1e3;
        C = 4.67
        %field = Monte.MW_power/(4*pi*(3e10/Monte.MW_Freq)*epsi*Volume)
        SIGMA0 = 1e-21
         dielectric_relax = 1e-8;
        
    end
    
    methods
      
    end
end

