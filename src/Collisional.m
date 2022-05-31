classdef Collisional
    %Collisional Class defining the reaction parameters of the system.
    
    properties
        PAPRES = 0;
        SIG1   = 0;
        WELL1  = 0;
        AMASS1 = 0;
        CVPI   = 0;
        CVPS   = 0;

        % 0.1	3.655	178.9	84.	5.00	0.00		! Krypton

        BGP     = 0;
        SIG2    = 0;
        WELL2   = 0;
        AMASS2  = 0;
        CVI     = 0;
        CVS     = 0;

        XSECT = 0;
        WELL  = 0;
        RMASS = 0;
        CV = 0;
        
        
        Coll_Step = zeros(2,10000);
        Cnormal = zeros(2,10000);
        
        Rate_1 = 10;
        Rate_2 = 10
        RRKM = zeros(1,10000);
    end
    
    
end

