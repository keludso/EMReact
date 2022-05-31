%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Monte Carlo Simulation of Master equation to include microwave effects
% on the rate of the reaction.
%
% Author - Kelvin Dsouza, Dr Daryoosh Vashaee 
%
% Comment - The code is is converted from the MULTIWELL code in fortran
%  which results in calculating the unimolecular dissociation by Dr Barker
%
% The code begins with setting the parameter s for the reaction to occur. 
% 1. Number of trials - Number of trajectories in the path of the reaction 
% 2. Max Energy - Energy limit for simulation default- 100000
% 3. Grain Size - Energy grain size that would be used in the simulation. 
%               smaller the grain size better resolution 
% 4. Time Limit - The max time for simulation 
% 5. Microwave ON/OFF - Choosing the microwave settings 
% 6. Microwave Power - The amount of Fluence deposited on the sample during
%                       simulation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; 

%% Input Configurations dialog box

prompt = {'File-name', 'Number of Trials','MAX Energy (cm-1)','Grid Size (cm-1):','Time Limit', ...
    'Microwave ON/OFF (1 or 0)', 'Microwave Power (J/cm2/s)', 'Starting Energy (cm-1)','Temperature (K)' ,'Frame size','Vibration data file (.txt)','Reaction Parameters (.txt)'};
 title = 'Energy Distribution parameters';
dims = [1 60];
definput = {'name_me.txt', '3000','50000','5', '5e-3', '1','2e3', '1000','300', '100','C6H6_vib.txt','reaction_parameters.txt' };
answer = inputdlg(prompt,title,dims,definput);

% Maximum energy and step size 

%% Global Variables 
% Assigning the Monte Carlo Class

cd src
Monte = MonteCarlo;
Microwave = Micro_wave;
Monte.answer = answer;
Monte.Emax = str2double(answer(3));
Monte.Step = str2double(answer(4));
Monte.Nmax = Monte.Emax/Monte.Step;  % Maximum energy = Nmax*Step
Energy = 1:Monte.Step:Monte.Emax;
Monte.Time_res = str2double(answer(10));

%% SET NUMBER OF TRIALS 
%
%  Nstart = NO. OF TRAJECTORIES  
%                         	     

Monte.Nstart = str2double(answer(2)); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%  Calculating the density of states   %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[density, sum] = calculate_density(Monte.Nmax,Monte.Step, answer{11});

load DMSO_dens.mat % I have used MULTIWELL to calculate the density of states
Monte.density_states = DMSO_dens;  
%Monte.ALNDEN = log(density);   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SET THE EBEGIN, FLAGS
% EBEGIN  - Starting energy point for simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Monte.Ebegin   = str2double(answer(8));

%% SETTING TEMPERATURE 
%   
%

Monte.Temp  = str2double(answer(9));
Monte.Temp_init   = Monte.Temp; % Temperature 


%% Energy Transfer Parameters 
%
%	Energy Transfer parameters
% DC Stores the coefficients to be used in the model

Monte.ITYPE   = [1,1];
Monte.DC(1,:) = [10, 0.017, -1.18e-07, 1.5e-3,10000,0.3,0,0];
Monte.DC(2,:) = [0.1, 0.00921, -1.728e-07, 0.1,100,0,0,0];
%Monte.DC(2,:) = [0.1, 0.00821, -0.838e-07, 0.1,100,0,0,0]; %worked well
%Monte.DC(2,:) = [10, 0.017, -1.18e-07, 1.5e-3,10000,0.3,0,0];
%Monte.DC(2,:) = [15.4, 0.00521, -0.738e-07, 0,0,0,0,0];

%% Field Parameters

%   Microwave parameters
%
%  
%

Microwave.MW_on      = str2double(answer(6));   %MW ON/OFF
Microwave.MW_power   = str2double(answer(7));
Microwave.MW_Freq    = 0.08; %2.4 GHz in cm-1

%%     TLIM IS THE TIME LIMIT OF THE CALCULATION
%8

Monte.TLIM   = str2double(answer(5));
Microwave.SIGMA0   = 1e-21;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulating the Monte Carlo Code 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Monte = Monte_Simulate_parfor(Monte,Microwave,answer{12}, answer{1});

% Processing the data 
% cd src 
% Monte = post_proc(Monte);

%% Writing to file 

write_output(Monte, Microwave);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
