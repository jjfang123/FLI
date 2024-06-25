function SelInd  = MRA(PopObj,IMatrix,W,N)
% calculate the size of individuals in the first nondominant level
IrFitness                         = min(IMatrix,[],2);
Level1Index                       = find(IrFitness>=0);
Len_Level1                        = length(Level1Index);

if Len_Level1<=N
    [~,SortIndex]                 = sort(-IrFitness);
    SelInd                        = SortIndex(1:N);
else
    % only focus on the solutions in the first level
    SelInd                        = Level1Index;
    PopObj                        = PopObj(Level1Index,:);
    IMatrix                       = IMatrix(Level1Index,Level1Index);
    
    [Num,M]                       = size(PopObj);
    NormW                         = W./repmat(sqrt(sum(W.^2,2)),1,M);
    NormPopObj                    = PopObj./repmat(sqrt(sum(PopObj.^2,2)),1,M);
    [~,ZoneIndex]                 = max(NormPopObj * NormW',[],2);
    Num_W                         = size(W,1);
    ZoneDensity                   = zeros(1,Num_W);
    zone.index                    = [];
    Zone                          = repmat(zone,1,Num_W);
    for j = 1:Num
        Zj                        = ZoneIndex(j);
        Zone(Zj).index            = [Zone(Zj).index,j];
        ZoneDensity(Zj)           = ZoneDensity(Zj) + 1;
    end
    
    [NDensity,SortIndex]          = sort(-ZoneDensity);
    Density                       = abs(NDensity);
    [Values,Neightboor]           = min(IMatrix,[],2);
    
    i                             = 1;
    RemoveIndex                   = [];
    while Num>N
        DelNum                    = min(Num - N,Density(i)-1);
        CandidateIndex            = Zone(SortIndex(i)).index;
        
        for j = 1:DelNum
            [~,Del_Ind]           = min(Values(CandidateIndex));
            Del_Ind               = CandidateIndex(Del_Ind);
            
            RemoveIndex           = [RemoveIndex,Del_Ind];
            IMatrix(Del_Ind,:)    = Inf;
            IMatrix(:,Del_Ind)    = Inf;
            Need_Updata           = find(Neightboor==Del_Ind);
            L_Need                = length(Need_Updata);
            if L_Need>0
                [Values(Need_Updata),Neightboor(Need_Updata)]=min(IMatrix(Need_Updata,:),[],2);
            end
            Values(Del_Ind)       = Inf;
        end
        Num                       = Num - DelNum;
        i                         = i+1;
    end
    
    SelInd(RemoveIndex)           = [];
end

end