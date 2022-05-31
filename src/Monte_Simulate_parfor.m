function Monte =  Monte_Simulate_parfor(Monte, Microwave, r_file, comment)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% MonteCarlo Simulation 
% Simulation of MonteCarlo code 
% 
% Finding starting energy in the distribution %%
% Runs the gillespie algorithm to calculate the energy distribution based on
% reaction rate 


%% Assigning variable 

% Collisonal Parameters 

Colli = Collisional;
Colli = read_reaction_param(Monte,Colli, r_file);
[Colli.Rate_1, Colli.Rate_2, Colli.RRKM] = calculate_rate(Monte,Colli);

% Normalization functions 
Monte.ICOLL = 1;
Colli.Coll_Step = zeros(2,10000);
Colli.Cnormal = zeros(2,10000);
Colli = Enormal(Monte,Colli);
Monte.ICOLL = 2;
Colli = Enormal(Monte,Colli);

TLIM = Monte.TLIM;
Time_res = 101;

% Initializing arrays 
EDIST = zeros(Time_res+1,Monte.Nmax);
Photon_absorbed = zeros(1,Time_res+1);
Photon_emitted = zeros(1,Time_res+1);
NN = zeros(Time_res+1,1);
NTRAJ_t =  cell(Monte.Nstart ,1);
Coll_time_t = cell(Monte.Nstart ,1);
absorbed_avg =zeros(1,Time_res);
emitted_avg =zeros(1,Time_res);
Temp_change = zeros(1,100);
NTRAJ = 0   ;
Coll_time = 0;
t = 1;
T_finish = 0;
Temp_change = zeros(Time_res,Monte.Nmax);
inital_Temp = Monte.Temp;
%% Running the loop

for j = 1:Monte.Nstart
   
    % Finding the starting energy 
      E = Monte.Ebegin;
      NET_A =0;
      NET_B =0;
      status = 0;
      T=0;
      Monte.Temp = inital_Temp;
      Monte.Temp_init = Monte.Temp;
      [Colli.Rate_1, Colli.Rate_2, Colli.RRKM] = calculate_rate(Monte,Colli);
      while (status == 0)
        
          EE = E ;      % Energy for bookkeeping
          TT = T;       % Time for Book keeping
          
          
          if T >= TLIM % If time has excedded the limit next itteration 
            T_finish = T_finish + 1;
            break;
          end     
          
          
          % Calculate the time step and probabilities
          
          
          
          [TAU, P] = Stochastic_calc(Monte,Microwave, Colli, EE);      
  
          if (TAU+ TT) > TLIM || isnan(TAU)
                TAU = 1.01*TLIM; 
          end
          
         % Book keeping function
        
         T_start = round(Time_res*TT / TLIM) + 2; %    DIVIDE TIME INTO 100 INTERVALS FROM TIME = 0.0 TO TLIM
         T_start(T_start>Time_res+1) = Time_res+1;
         AARGU = (TT + TAU)/TLIM;
         T_finish = round(Time_res.*AARGU) + 2;
         T_finish(T_finish>Time_res+1) = Time_res+1;
         
         K = 1 + round(EE/Monte.Step);
     
        if K > length(Monte.density_states)
            K = length(Monte.density_states);
        end
        
        if TT == 0   %	AT T = 0 
            EDIST(1,K) = EDIST(1,K) + 1;
            NN(1) = NN(1) + 1;
        end

        NN(T_start:T_finish) = NN(T_start:T_finish) + 1;
        EDIST(T_start:T_finish,K) = EDIST(T_start:T_finish,K) + 1;
        Photon_absorbed(T_start:T_finish) = Photon_absorbed(T_start:T_finish) + NET_A;
        Photon_emitted(T_start:T_finish) = Photon_emitted(T_start:T_finish) - NET_B;
        %Temp_change(T_start:T_finish)
        
        % Adjusting for Temperature 
        if Monte.Temp_init + 1 < Monte.Temp            
            Monte.Temp_init = Monte.Temp;
            T_step = round(100*TT/TLIM)+1;
            Temp_change(T_step,j) = Monte.Temp;
            [Colli.Rate_1, Colli.Rate_2, Colli.RRKM] = calculate_rate(Monte,Colli);
        end
        
        
        if TAU == 0
            status = 0;
        elseif T > TLIM
            status = 1;
        else
            T = T+ TAU;
            [status,E, NET_A, NET_B, Monte] = monte_calc(Monte, Colli,EE, P,TAU, Microwave);
                            
            %Monte.Temp = Temp + MW_temp;
        end  
        
        if E <= 0
            E = Monte.Step;
        end
        if K < 10
            pause(1);
        end
      end
      
      if T < TLIM
          NTRAJ(j) = 1;
          Coll_time(j) = TT;
      end
      NN(NN==0) = 1;

  j   
end   
 

% Calculating Trajectories 
sum_r = 0;
Coll = zeros(1,101);
for tm=1:Time_res+1
        for k = 1:length(Coll_time)
            if tm == round(100*Coll_time(k)/TLIM)
                sum_r = sum_r +1;
            end
        end
    Coll(tm) = sum_r;
    
end

Energy = 1:Monte.Step:Monte.Emax; %max
Time = 0:TLIM/Time_res:TLIM;
%Edistribution = Monte.Nstart*avg_energy/j;

% Saving to data structure 

Monte.Photon_absorbed = Photon_absorbed;
Monte.Photon_emitted = Photon_emitted;
Monte.Edistribution = (EDIST./NN)*Monte.Nstart;
Monte.Traj_reacted = Coll;
Monte.Time = Time;
Monte.Energy = Energy;
Monte.Temp_change = Temp_change;
Monte.Coll = Coll;


return
end

