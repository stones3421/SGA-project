%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The script implements the template for the comparison algorithm%%%%%%%%%%%
%%%% %%%%%%%%Author by: Ai-Qing Tian%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% %%%%%%%%All Rights Reserved.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Unauthorized copying of this file, via any medium is strictly prohibited%%%%
%%%%Copyright 2022-07-16 username stones12138@163.com%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all 
clc

SearchAgents_no=30; % Number of search agents

Function_name='F8'; % Name of the test function that can be from F1 to F23 (Table 1,2,3 in the paper)

Max_iteration=1000; % Maximum numbef of iterations

% Load details of the selected benchmark function
[lb,ub,dim,fobj]=Get_Functions_details(Function_name);

[Best_score_sga,Best_pos_sga,cg_curve_sga]=SGA(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);% SGA算法


figure('Position',[284   214   660   290])
%Draw search space
subplot(1,2,1);
func_plot(Function_name);
title('Test function')
xlabel('x_1');
ylabel('x_2');
zlabel([Function_name,'( x_1 , x_2 )'])
grid off

%Draw objective space
subplot(1,2,2);
semilogy(cg_curve_sga,'Color','b')
title('Convergence curve')
xlabel('Iteration');
ylabel('Best flame (score) obtained so far');

axis tight
grid off
box on
legend('SGA')

display(['The best solution obtained by SCA is : ', num2str(Best_pos_sga)]);
display(['The best optimal value of the objective funciton found by SCA is : ', num2str(Best_score_sga)]);