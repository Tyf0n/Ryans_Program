clear;

% This is necessary for Tower usage
cd '/home/tower/rpeng/Newer_Program_CPLEX';
addpath('/opt/ilog/cplex/matlab/');

% Opening the parallel session
if matlabpool('size')== 0
    matlabpool open
end

tic % Start timer

B = 10000000;
epsilon = 0.001;

% Main variables and parameters
s0 = 0;  % Initial state
t  = 1;  % Starting time point
T  = 5 * 365; % Number of days and depth of simulation
threshold = 10; % Threshold value
eta  = 1; % Discount factor for the risk averse part of the measure
Nsim = 40;  % This is needed just for the moment while we deal with the y's
demand = 10000; % This is just used for the moment to replace the Lpathselectedhour with some fixed demand

% folder where the output will be saved
parentFolder = 'OutputData';
folder =  sprintf('T%u_th%g_eta%g_Nsim%u',T,threshold,eta,Nsim);
pathToFolder = sprintf('OutputData/T%u_th%g_eta%g_Nsim%u',T,threshold,eta,Nsim);
[status,message,messageid] = mkdir(parentFolder, folder);

% diary(sprintf('%s/log.txt',pathToFolder));

% Write the main variables and parameters to the log
s0
t
T
threshold
eta
Nsim
demand

% Load data values
load 'DataGeneration/Data_variables.mat'
% load 'Lpathselectedhour.mat';

YPathselectedhour = zeros(T,Nsim);
Lpathselectedhour = demand * ones(T,Nsim);

outputdata =  zeros(Nsim,2);

for trend = 1:8
    
    switch trend
        case 1
            Ppathselectedhour = Pgenerated_uplin;
            Fselectedhour     = Fgenerated_uplin;
        case 2
            Ppathselectedhour = Pgenerated_downlin;
            Fselectedhour     = Fgenerated_downlin;
        case 3
            Ppathselectedhour = Pgenerated_upquad;
            Fselectedhour     = Fgenerated_upquad;
        case 4
            Ppathselectedhour = Pgenerated_downquad;
            Fselectedhour     = Fgenerated_downquad;
        case 5
            Ppathselectedhour = Pgenerated_sine;
            Fselectedhour     = Fgenerated_sine;
        case 6
            Ppathselectedhour = Pgenerated_uposc;
            Fselectedhour     = Fgenerated_uposc;
        case 7
            Ppathselectedhour = Pgenerated_downosc;
            Fselectedhour     = Fgenerated_downosc;
        case 8
            Ppathselectedhour = Pgenerated_rotatelin;
            Fselectedhour     = Fgenerated_rotatelin;
    end
    
    for type = 1:4
        
        parfor i = 1:Nsim
            [x{i},V{i},Vbar{i},fval{i},exitflag{i},outputdatacell{i}] = ...
                LPSolver(s0,t,T,threshold,eta,i,type,Ppathselectedhour,Lpathselectedhour,Fselectedhour,YPathselectedhour);
        end
        
        % Rearrange output data into structured form to allow for post-processing
        for j = 1:Nsim
            outputdata(j,1) = outputdatacell{j}(1);
            outputdata(j,2) = outputdatacell{j}(2);
        end
        
        % Get rid of NaNs in the output data
        outputdata(any(isnan(outputdata),2),:) = [];

        % Delete rows that have maxed out x values and replace with a single row at
        % the end with the average value.
        maxedrows = outputdata(outputdata(:,1) == B,:);
        mrlength = length(maxedrows);

        if mrlength ~= 0
        
            mrtotal = sum(maxedrows);
            mraverage = mrtotal / mrlength;

            outputdata((outputdata(:,1) >= B-epsilon & outputdata(:,1) <= B+epsilon),:) = [];

            outputlength = length(outputdata);
            outputdata(outputlength+1,:) = mraverage;

        end
        
        % Create new variables for name-stamping
        eval(sprintf('outputdata_trend%d_type%d = outputdata;',trend,type));
        eval(sprintf('x_trend%d_type%d = x;',trend,type));
        eval(sprintf('fval_trend%d_type%d = fval;',trend,type));
        eval(sprintf('V_trend%d_type%d = V;',trend,type));
        eval(sprintf('Vbar_trend%d_type%d = Vbar;',trend,type));
    end
end

elapsedtime = toc; % Stop timer

% Save the output
save(sprintf('%s/Output_Variables.mat',pathToFolder), 'outputdata_*','x_*','fval_*','V_*','Vbar_*','elapsedtime');

% diary off;

matlabpool close % Closing parallel session

% Testing the Github branching system!
