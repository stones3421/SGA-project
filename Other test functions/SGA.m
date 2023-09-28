function [gBestScore,gBest,cg] = SGA6(Np, Max_iter, Xmin, Xmax, D, fobj)

Pos = initialization(Np,D,Xmax,Xmin);
Vel = zeros(Np,D);

gBest = zeros(1,D);
gBestScore = inf;
cg = zeros(1,Max_iter);

for i = 1:Np
    fitness(i) = fobj(Pos(i,:));
end

[gBestScore, index] = min(fitness);
gBest = Pos(index,:);

cg(1) = gBestScore;

t = 2;

while t < Max_iter
    coe = (4*(t/Max_iter))/exp(4*(t/Max_iter));
    fi = rand()*2*pi;
    for i = 1:Np
        acc = ((gBest - Pos(i,:)) - 1.29*Vel(i,:).^2*sin(fi))*10^-2;
        Vel(i,:) = coe*Vel(i,:) + acc;
    end

   [~,index] = sort(fitness);
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
        aa(i,:) = Pos(i,:)*fitness(i);
        bb(i) = Np*fitness(i);
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
                Flag4ub=Pos(i,:)>Xmax;
                Flag4lb=Pos(i,:)<Xmin;
                Pos(i,:)=(Pos(i,:).*(~(Flag4ub+Flag4lb)))+Xmax.*Flag4ub+Xmin.*Flag4lb;
            end
    else
        if rand>0.5
            for i = 1:Np
                Pos(i,:) = Pos(i,:) + (Pos(i,:) - gBest)*rand;
                Flag4ub=Pos(i,:)>Xmax;
                Flag4lb=Pos(i,:)<Xmin;
                Pos(i,:)=(Pos(i,:).*(~(Flag4ub+Flag4lb)))+Xmax.*Flag4ub+Xmin.*Flag4lb;
            end
        else
            for i = 1:Np
                Pos(i,:) = gBest + (Pos(i,:) - gBest).*rand.*Brownian(D);
                Flag4ub=Pos(i,:)>Xmax;
                Flag4lb=Pos(i,:)<Xmin;
                Pos(i,:)=(Pos(i,:).*(~(Flag4ub+Flag4lb)))+Xmax.*Flag4ub+Xmin.*Flag4lb;
            end
        end
    end        

    for i = 1:Np
        fitness(i) = fobj(Pos(i,:));
    end

    for i=1:Np
        if gBestScore > fitness(i)
            gBest = Pos(i,:);
            gBestScore = fitness(i);
        end
    end
    t = t +1;
    cg (t) = gBestScore;
end
end


function o = Brownian(D)
    T = 1;
    r = T/D;
    dw = sqrt(r)*randn(1,D);
    o = cumsum(dw);
end
