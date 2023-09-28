function [gBest,gBestScore,cg,FES] = SGA(Np, D, Pos, Vel, Xmax, Xmin, fhd, Max_Fes,varargin)
%% Snow Geese Algorithm
% This method implements the snow goose algorithm
% input：Np-->Population size
%               D  -->Dimension of the problem
%               Pos-->Location of initial stock
%               fhd-->Selection of issues
%               Max_Fes-->Maximum number of evaluations
%               varargin-->Selection of functions
%output：gBest -->optimal candidate solution
%           gBestScore-->optimum value
%           cg              -->Iterative convergence curves
%           FES           -->Number of evaluations
%% This algorithm implements the intelligent behaviour of the process of 
% snow geese flying over long distances and contains two main aspects

%% Algorithm Start
%Initial evaluation of algorithms
if size(Xmax,1) == 1
    Xmax = repmat(Xmax,1,D);
    Xmin = repmat(Xmin,1,D);
end
e=feval(fhd,Pos',varargin{:});

FES = Np;
[gBestScore,index] = min(e);
gBest = Pos(index,:);
iter = 1;
Max_Iter = Max_Fes/Np;
cg = zeros(1,Max_Iter);
cg(1) = gBestScore;

%enter main loop
while FES < Max_Fes
    coe = (4*(iter/Max_Iter))/exp(4*(iter/Max_Iter));
%     coe = rand;
%          coe = iter/Max_Iter;
%     fi = (iter/Max_Iter)*2*pi;s
    fi = rand()*2*pi;
    for i = 1:Np
        acc = ((gBest - Pos(i,:)) - 1.29*Vel(i,:).^2*sin(fi))*10^-2;
        Vel(i,:) = coe*Vel(i,:) + acc;
    end

   [~,index] = sort(e);
   for i = 1: Np
       New_Pos(i,:) = Pos(index(i),:);
       New_Vel(i,:) = Vel(index(i),:);
   end
   Pos = New_Pos;
   Vel = New_Vel;
    a = 4*rand() - 2;
    b = 3*rand() -1.5;
    c = 2*rand() - 1;

    for i =1:Np
        aa(i,:) = Pos(i,:)*e(i);
        bb(i) = Np*e(i);
    end
    Xc = sum(aa)/sum(bb);
    Pos = Pos + Vel;
    if fi < pi
            for i = 1:Np
                if i<=1/5*Np
                    Pos(i,:) = Pos(i,:)  + a*(gBest - Pos(i,:)) + Vel(i,:);    
                elseif  1/5*Np< i && i < 4/5*Np
                    Pos(i,:) = Pos(i,:)  + a*(gBest - Pos(i,:)) + b*(Xc -Pos(i,:)) - c*(Pos(Np,:) + Pos(i,:))  + Vel(i,:);  
                else
                    Pos(i,:) = Pos(i,:)  + a*(gBest - Pos(i,:)) + b*(Xc -Pos(i,:)) + Vel(i,:);  
                end
            end
    else
        if rand>0.5
            for i = 1:Np
                Pos(i,:) = Pos(i,:) + (Pos(i,:) - gBest)*rand;
            end
        else
            for i = 1:Np
                Pos(i,:) = gBest + (Pos(i,:) - gBest).*rand.*Brownian(D);
            end
        end
    end


    for i =1:Np
        for j = 1:D
            if Pos(i,j) > Xmax(j)
                Pos(i,j) = Xmax(j);
            end
            if Pos(i,j) < Xmin(j)
                Pos(i,j) = Xmin(j);
            end
        end
    end

    e = feval(fhd,Pos',varargin{:});
    FES = FES +Np;

    for i=1:Np
        if gBestScore > e(i)
            gBest = Pos(i,:);
            gBestScore = e(i);
        end
    end
    iter = iter +1;
    cg (iter) = gBestScore;
end
end


function o = Brownian(D)
    T = 1;
    r = T/D;
    dw = sqrt(r)*randn(1,D);
    o = cumsum(dw);
end
