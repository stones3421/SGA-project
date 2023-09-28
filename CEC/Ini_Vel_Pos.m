function [Pos,Vel] = Ini_Vel_Pos(Xmax, Xmin, Np, D)
%This function implements the process of initialising the algorithm.
    rand('seed',sum(100*clock));
    Pos = Xmin + (Xmax - Xmin).*lhsdesign(Np,D);%
    Vel = Xmin/10 + (Xmax - Xmin)/10.*lhsdesign(Np,D);%
end