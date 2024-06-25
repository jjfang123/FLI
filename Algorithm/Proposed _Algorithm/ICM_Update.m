function [Population,Archive] = ICM_Update(MaxPop,N,W,Zmin,Fmin)

CV                   = sum(max(0,MaxPop.cons),2);
PopObj               = MaxPop.objs;
[Num,M]              = size(PopObj);
% shift the objective space to R+
PopObj               = PopObj - repmat(Zmin,Num,1) + 1e-6;

% calculate the indicator matrix
IMatrix              = ones(Num,Num);

for i = 1:1:Num
    Fi               = PopObj(i,:);
    Ci               = CV(i);
    
    if Ci == 0 %%%%% Xi is feasible
        Ir               = log(repmat(Fi,Num,1)./PopObj);
        MaxIr            = max(Ir,[],2);
        MinIr            = min(Ir,[],2);
        CVA              = MaxIr;
        DomInds          = find(MaxIr<=0);
        CVA(DomInds)     = MinIr(DomInds);
    else  %%%%% Xi is an infeasible solution
        IC               = log(repmat(Ci+1e-6,Num,1)./(CV+1e-6));
        Fi               = repmat(Fi,Num,1);
        CVA              = max(log(max(max(Fi,PopObj)./min(Fi,PopObj),[],2)),IC);
    end
    
    IMatrix(:,i)     = CVA;
    IMatrix(i,i)     = Inf;
end

FeasibleInd                   = find(CV==0);
Len_F                         = length(FeasibleInd);
if Len_F<=N
    [~,CV_SortInd]            = sort(CV);
    Archive                   = MaxPop(CV_SortInd(1:N));
else
    FPopObj                   = PopObj(FeasibleInd,:) + repmat(Zmin,Len_F,1) - repmat(Fmin,Len_F,1);
    SelInd                    = MPRPD(FPopObj,IMatrix(FeasibleInd,FeasibleInd),N,1);
    Archive                   = MaxPop(FeasibleInd(SelInd));
end

% %%%%%%%%%%%%%%%   一般的方法
% SelInd                        = MPRPD(PopObj,IMatrix,N,2);
% Population                    = MaxPop(SelInd);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%  a mechanism based on resource assign
SelInd                        = MRA(PopObj,IMatrix,W,N);
Population                    = MaxPop(SelInd);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end