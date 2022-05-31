function Monte = post_proc(Monte)

% This function calculates the variance of the plot by fitting a gaussian 
% curve to the input. 

% Finding the mean of the energy distribution by finding the maximum point 

Energy_distr = Monte.Edistribution;
Energy = Monte.Energy;
temp_count = 1;
RMSE_O = 1;
% Looping through all stages of time
c=1;

for i =80:101
%% Fit:

edi = log(smoothdata(Energy_distr(i,20:4000),'gaussian','SmoothingFactor',0.05));
[xData, yData] = prepareCurveData( Energy(20:4000), edi );

% Set up fittype and options.
ft = fittype( 'a*log(x) -1.5*log(c) -x/c + d', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 -Inf];
opts.StartPoint = [0.0288029447581492 0.525109224225604 0.0982937021602779];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

Monte.kT(c) = fitresult.c;
c = c + 1;
end

% figure;
% % 
% y1 = exp(fitresult.a.*log(Energy)).*exp(-1.5*log(fitresult.c)).*exp(-Energy./fitresult.c).*exp(fitresult.d);
% 
% plot(Energy,y1);
% hold on
% % %plot(Energy,y2);
% % %hold on
% plot(Energy, smooth(smooth(Edistribution(100,:))));

% % figure;
% plot(Monte.Variance_temp);
% mean_T = mean(Monte.Variance_temp);
% hold on;
% plot(ones(length(Monte.Variance_temp),1)*mean_T)

% Fitting through log 

%log(y) = log(a*E) - E/kT


    %   T =        1623  (1407, 1840)
     %  a =     0.02543  (0.01371, 0.03715)
     
     
% distri = 0.03543.*Energy.*exp(-Energy./1623);
% figure;
% plot(Energy, distri);
% hold on;
% plot(Energy,Monte.Edistribution(101,:));





end

