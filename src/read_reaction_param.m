function [Collisional] = read_reaction_param(Monte,Collisional,file_name)



%% SETTING COLLISION PARAMETERS 


% Details :
% The file should be in a txt file format and the following structure 
% 
% Line 1 :    "<Title>"
% Line 2 : 'Empty' Description for the user  
% Line 3 : <Number of vibrational mode>  <Number of Rotational Modes>
% Line 4 : <Mode = 'vib'/'rot'  <Frequency/ Ie>
% Line 5 : <Mode = 'vib'/'rot'  <Frequency/ Ie>
%

%%%%%%%% Constant Defination %%%%%%%%%%%%%%%%%%%

cd ..
cd input_files

rk_file = fopen(file_name);
i = 1;
tline = fgetl(rk_file);
title = tline;
comment = fgetl(rk_file);
tline = fgetl(rk_file);
i = 1;
while i<=16
    
   tline = fgetl(rk_file); 
   Str_file= tline;
   splt = strsplit(Str_file,' ');
   LJ_Param(i) = str2double(splt{2});
   i = i+1;
   
end 

fclose(rk_file);
cd ..
cd src

%% Assigning the values 
Collisional.PAPRES = LJ_Param(1);
Collisional.SIG1   = LJ_Param(2);
Collisional.WELL1  = LJ_Param(3);
Collisional.AMASS1 = LJ_Param(4);
Collisional.CVPI   = LJ_Param(5);
Collisional.CVPS   = LJ_Param(6);

% 0.1	3.655	178.9	84.	5.00	0.00		! Krypton

Collisional.BGP     = LJ_Param(7);
Collisional.SIG2    = LJ_Param(8);
Collisional.WELL2   = LJ_Param(9);
Collisional.AMASS2  = LJ_Param(10);
Collisional.CVI     = LJ_Param(11);
Collisional.CVS     = LJ_Param(12);

Collisional.XSECT = 0.5*(Collisional.SIG1 + Collisional.SIG2);
Collisional.WELL = sqrt(Collisional.WELL1*Collisional.WELL2);
Collisional.RMASS = Collisional.AMASS1*Collisional.AMASS2/(Collisional.AMASS1 + Collisional.AMASS2);
Collisional.BGP = Collisional.BGP*9.66e18 /Monte.Temp;
Collisional.PAPRES = Collisional.PAPRES*9.66e18/Monte.Temp;
Collisional.CV = (Collisional.BGP*(Collisional.CVI+(Monte.Temp - 300)*Collisional.CVS) + Collisional.PAPRES*(Collisional.CVI + (Monte.Temp - 300)*Collisional.CVPS))*5.808e-25;



end