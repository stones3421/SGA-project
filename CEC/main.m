%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The script implements the template for the comparison algorithm%%%%%%%%%%%
%%%% %%%%%%%%Author by: Ai-Qing Tian%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% %%%%%%%%All Rights Reserved.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Unauthorized copying of this file, via any medium is strictly prohibited%%%%
%%%%Copyright 2022-07-16 username stones12138@163.com%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all 
clc
warning off

% Format of results
format long;

% Initialising test function properties
fhd = str2func('cec14_func');
Fnum = 30;

D = 30; % Problem testing dimensions
Np = 120; % Population size
Max_Fes = Np*10000; % Set the maximum number of evaluations according to CEC2014 test document requirements
Xmax = 100; % Searching the Upper boundary
Xmin = -100; %  Searching the Lower boundary

for f = 1:Fnum
    [Pos,Vel] = Ini_Vel_Pos(Xmax,Xmin,Np,D); % Initial position and speed
    [gBest_SGA,gBestScore_SGA,cg_SGA,FES_SGA] = SGA(Np, D, Pos, Vel, Xmax, Xmin, fhd, Max_Fes,f);
    figure;
    plot(cg_SGA,'Color','b')
    title('Convergence curve of function', num2str(f));
    xlabel('Iteration');
    ylabel('Best flame (score) obtained so far');
end

