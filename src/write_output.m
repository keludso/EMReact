function write_output(Monte,Microwave)

%%Writing the values from the class to a file for storage

cd ..
cd data_file

output_file = fopen(Monte.answer{1},'w');
fprintf(output_file,'########################################################\n');
fprintf(output_file ,'##### Internal Energy Evaluation Tool ##################\n');
fprintf(output_file,'########################################################\n\n\n');

fprintf(output_file,'    This program is a version of Dr Barkers Multiwell Code\n');
fprintf(output_file,'    which has been written in MATLAB\n');
fprintf(output_file,'    Please cite ....\n\n');
fprintf(output_file,'    This program is a version of Dr Barkers Multiwell Code\n');

fprintf(output_file,'--- Input Specifications ----\n\n');
fprintf(output_file,'Number of Trajectories# %d\n', Monte.Nstart);
fprintf(output_file,'Time Limit Set(seconds)# %f\n', Monte.TLIM);
fprintf(output_file,'Time Resolution # %d seconds\n', Monte.Time_res);
fprintf(output_file,'Emax # %d \n', Monte.Emax);
fprintf(output_file,'Energy Step # %d\n', Monte.Step);
fprintf(output_file,'Microwave - #%d\n', Microwave.MW_on );
fprintf(output_file,'Microwave Power #%d\n', Microwave.MW_power);
fprintf(output_file,'Microwave Frequency(cm-1) #%f\n', Microwave.MW_Freq );
% Number of Trials
% Time ,
% MW Switch 
% MW power 

fprintf(output_file,'--- Final Results ----\n');
% Trajectories Reacted 
% Photon Absorbed 
% Photon Emitted 
% Reaction


% Density of States 


%Energy distribution 

fprintf(output_file,'Energy distribution \n');


% Writing the energy distribution

for edd= 1:Monte.Time_res+1
    
    A = Monte.Edistribution(edd,:);
    fprintf(output_file,'%d ', A);
    fprintf(output_file,'\n ');
end

% Photon absorbed 

fprintf(output_file,'\n\n ');
fprintf(output_file,'Photon_absorbed\n ');
A = Monte.Photon_absorbed;
fprintf(output_file,'%d ', A);
fprintf(output_file,'\n\n');
fprintf(output_file,'emitted\n ');
A = Monte.Photon_emitted;
fprintf(output_file,'%d ', A);
fprintf(output_file,'\n\n');

fprintf(output_file,'Reaction \n');
fprintf(output_file,'%d ', Monte.Coll);
fprintf(output_file,' \n');
fprintf(output_file,'Temperature change\n');
for edd= 1:Monte.Time_res+1
    
    A = Monte.Temp_change(edd,:);
    fprintf(output_file,'%d ', A);
    fprintf(output_file,'\n ');
end

fclose(output_file);

cd ..
cd src




end

