classdef MonteCarlo
    % A library for constants and functions 
    
    properties
        
        %% Input Parameters 
        answer = cell(12,1);
        
        Step =0;            % Grain Size
        TLIM = 0;
        Emax =0;
        Nmax =0;
        density_states = zeros(1,4000);
        h = 1.9863e-23;
        Temp, TempI = 0;
        Ebegin = 0;
        Nstart = 0;
        DC = zeros(2,8);
        ITYPE   = [1,1];
        MW_power = 0;
        SIGMA0 = 0;
        Time_res =0;
        MW_on = 0;
        MW_Freq = 0;
        MW_absorbed = 0;
        Energy = 0;
        ICOLL = 0;
        Temp_change = 0;
        
        % Output Parameters 
        Traj_reacted = 0;
        Edistribution= 0;
        Photon_absorbed = 0;
        Photon_emitted = 0;
        Time = 0;
        kT = 0;
        Temp_calc = 0;
        IFLAG =0;
        Temp_init = 0;
        Coll= 0;
        
    end
    
end

