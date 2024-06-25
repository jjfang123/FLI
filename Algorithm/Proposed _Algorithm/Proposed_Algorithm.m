function Proposed_Algorithm(Global)
% <algorithm> <A>
% Promising-region based EMO algorithm

%------------------------------- Reference --------------------------------
% Jiawei Yuan, Hai-Lin Liu,Fangqing Gu, Qingfu Zhang and Zhaoshui He,
% Investigating the properties of indicators and an evolutionary 
% many-objective algorithm based on a promising region. IEEE
% Transactions on Evolutionary Computation, 2020.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2020 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

%% Generate the random population
Population   = Global.Initialization();
Zmin         = min(Population.objs,[],1);
Fmin         = min(Population(all(Population.cons<=0,2)).objs,[],1);
Archive      = Population;
W            = UniformPoint(Global.N,Global.M);
Nr           = 1;
Ind          = calInd(Population,Global.gen,Global.maxgen);

%% Optimization
while Global.NotTermination(Archive)
    MixNum     = floor(Nr*Global.N);
    MatingPool = [Population(randsample(Global.N,MixNum)),Archive(randsample(Global.N,Global.N-MixNum))];
    if mod(Global.gen,100)==0
        Ind  = calInd(Population,Global.gen,Global.maxgen);
    end
%    MatingPool = Population(randsample(Global.N,Global.N));
 if rand<Ind
     if rand>0.5
        Offspring  = DE(MatingPool,MatingPool(randperm(Global.N)),Archive(randperm(Global.N)),{0.5,0.5,0.5,0.5});
     else
        Offspring  = DE(MatingPool,MatingPool(randperm(Global.N)),Archive(randperm(Global.N)),{1,1,1,1});
     end
 else

      %  Offspring  = DE(Population,Population(randperm(Global.N)),Population(randperm(Global.N)));
      %   Offspring  = DE_cr(Population,Population(randperm(Global.N)),Population(randperm(Global.N)),Population(randperm(Global.N)),{1,0.5,1,20});
       Offspring  = DE_cr(MatingPool,MatingPool(randperm(Global.N)),MatingPool(randperm(Global.N)),MatingPool(randperm(Global.N)),{1,0.5,1,20});

end
%     
    Fmin       = min([Fmin;Offspring(all(Offspring.cons<=0,2)).objs],[],1);
    Zmin       = min([Zmin;Offspring.objs],[],1);
    [Population,Archive] = ICM_Update([Population,Offspring,Archive],Global.N,W,Zmin,Fmin);
    
    Nr         = max(1 - Global.gen/Global.maxgen,0.05);

end
end