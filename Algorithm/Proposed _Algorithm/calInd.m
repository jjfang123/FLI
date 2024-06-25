function Ind = calInd(Population,gen,maxgen)
%CALIND calculate the indicator according to the population and gen

N = 20; % number of reference points

PopObj          =  Population.objs;


%% Calculate the cons
Popcon               = sum(max(0,Population.cons),2);

%% select the next archive
FeasibleInd          = find(Popcon==0);
Len_F                = length(FeasibleInd);

P               =  Population(FeasibleInd);
PopObj          =  P.objs;
% A = size(unique(PopObj,'rows'),1);
% B = 1.1*N;
% if A<=B
%     N = floor(size(unique(PopObj,'rows'),1)/2);
%     Population = unique(Population);
% end


if Len_F < N
    [~,cons_index]   = sort(Popcon);
    Archive          = Population(cons_index(1:N));
    Rpop             = Population(cons_index(N+1:end));
else
    
%     Population       = Population(FeasibleInd);
% %    Fitness          = CalFitness(Population);
% %    Archive          = Environment(Population,Fitness,N);  
%     PopObj          =  Population.objs;
%     Zmin            =  min(PopObj,[],1);
%     NextIndex = MPRP(PopObj,Zmin,N);
%     Archive   = Population(NextIndex);
%     Rpop      = Population(setdiff([1:size(Population,2)]',NextIndex));
    
    
%%改动

    PopObj          =  Population.objs;
    Zmin            =  min(PopObj,[],1);
    NextIndex = MPRP(PopObj,Zmin,N);
    Archive   = Population(NextIndex);
    Rpop      = Population(setdiff([1:size(Population,2)]',NextIndex));
end

%% calculate the distance from each rest point to the population  
Adec = Archive.decs;
Rdec = Rpop.decs;
Dmatrix = zeros(size(Rdec,1),size(Adec,1));
for i = 1:size(Rdec,1)
    for j = 1:size(Adec,1)
        Dmatrix(i,j) = sum((Rdec(i,:)-Adec(j,:)).^2);
    end
end
Dis     = min(Dmatrix,[],2);

Zmin    = min(Rpop.objs,[],1);
Fit     = fitness(Rpop.objs,Zmin);
% [~,Fitrank] = sort(Fit);
% [~,Disrank] = sort(Dis);

R = corrcoef(Fit,Dis); % pearson corrcoef: [-1,1]
if size(R,1)==1
   B = NaN; 
end
R = R(1,2);

Ind = (0.5*R+0.5)^0.5;%平移到0-1区间
A = gen/maxgen;
Ind = Ind*A;

end

function NextIndex = MPRP(PopObj,Zmin,N)
[Num,M]                  = size(PopObj);
%% shift the objective space to R+
PopObj               = PopObj - repmat(Zmin,Num,1) + 1e-6;
%% calculate the indicator matrix
IMatrix              = ones(Num,Num);
for i = 1:1:Num
    Fi               = PopObj(i,:);
    Ir               = log(repmat(Fi,Num,1)./PopObj);
    MaxIr            = max(Ir,[],2);
    MinIr            = min(Ir,[],2);
    CVA              = MaxIr;
    DomInds          = find(MaxIr<=0);
    CVA(DomInds)     = MinIr(DomInds);
    IMatrix(:,i)     = CVA;
    IMatrix(i,i)     = Inf;
end

% calculate the size of individuals in the first nondominant level
IrFitness                = min(IMatrix,[],2);
Level1Index              = find(IrFitness>=0);
Len_Level1               = length(Level1Index);

if Len_Level1<=N
    [~,SortIndex]    = sort(-IrFitness);
    NextIndex        = SortIndex(1:N);
else
    % only focus on the solutions in the first level
    AllIndex         = Level1Index;
    PopObj           = PopObj(Level1Index,:);
    IMatrix          = IMatrix(Level1Index,Level1Index);
    
    %% select the valuable solutions in the current population
    MiddleIMatrix    = IMatrix;
    Ag1_n            = Len_Level1 - N;
    [Values,Neightboor] = min(MiddleIMatrix,[],2);
    BestInd          = 1:Len_Level1;
    Have_Delect      = zeros(1,Ag1_n);
    
    % mark the N solutions with the best fitness value by excluding the
    % worst individuals one by one
    for i=1:Ag1_n
        [~,Del_Ind]  = min(Values);
        Have_Delect(i) = Del_Ind;
        MiddleIMatrix(Del_Ind,:) = Inf;
        MiddleIMatrix(:,Del_Ind) = Inf;
        Need_Updata = find(Neightboor==Del_Ind);
        L_Need=length(Need_Updata);
        if L_Need>0
            [Values(Need_Updata),Neightboor(Need_Updata)]=min(MiddleIMatrix(Need_Updata,:),[],2);
        end
        Values(Del_Ind) = Inf;
    end
    BestInd(Have_Delect) = [];
    
    % determine the boundary of promising region
    Zmax                 = max(PopObj(BestInd,:),[],1);
    
    % remove the individuals outside the promising region
    OutIndex             = find(min(repmat(Zmax,Len_Level1,1) - PopObj,[],2) < 0);
    AllIndex(OutIndex)   = [];
    PopObj(OutIndex,:)   = [];
    IMatrix(OutIndex,:)  = [];
    IMatrix(:,OutIndex)  = [];
    Num                  = length(AllIndex);
    
    
    %% diversity maintance mechanism based on parallel distance
    % normalize the promising region
    PopObj               = PopObj./repmat(Zmax,Num,1);
    
    [Ir_Values,Ir_Neightboor] = min(IMatrix,[],2);
    
    DelectInd2 = [];
    
    DMatrix              = zeros(Num,Num);
    
    % calculate the parallel distance matrix
    for i = 1:Num
        Fi = PopObj(i,:);
        Fdelta = PopObj - repmat(Fi,Num,1);
        DMatrix(i,:) = sqrt(sum(Fdelta.^2,2) - (sum(Fdelta,2)).^2./M);
        DMatrix(i,i) = Inf;
    end
    
    [Dis_Values,Dis_Neightboor] = min(DMatrix,[],2);

    for l=1:(Num - N)
        [~,individual1]=min(Dis_Values);
        individual2=Dis_Neightboor(individual1);
        
        if Ir_Values(individual1)<Ir_Values(individual2)
            k = individual1;
        else
            k = individual2;
        end
        
        DelectInd2 = [DelectInd2,k];
        
        DMatrix(k,:) = Inf;
        DMatrix(:,k) = Inf;
        Need_Updata_Dis=find(Dis_Neightboor==k);
        L_Need_Dis=length(Need_Updata_Dis);
        if L_Need_Dis>0
            [Dis_Values(Need_Updata_Dis),Dis_Neightboor(Need_Updata_Dis)]=min(DMatrix(Need_Updata_Dis,:),[],2);
        end
        Dis_Values(DelectInd2)= Inf;

        IMatrix(k,:) = Inf;
        IMatrix(:,k) = Inf;
        Need_Updata_Ir=find(Ir_Neightboor==k);
        L_Need_Ir=length(Need_Updata_Ir);
        if L_Need_Ir>0
            [Ir_Values(Need_Updata_Ir),Ir_Neightboor(Need_Updata_Ir)]=min(IMatrix(Need_Updata_Ir,:),[],2);
        end
    end
    AllIndex(DelectInd2) = [];
    
    NextIndex = AllIndex;
end

end

function IrFitness = fitness(PopObj,Zmin)
[Num,M]                  = size(PopObj);
%% shift the objective space to R+
PopObj               = PopObj - repmat(Zmin,Num,1) + 1e-6;
%% calculate the indicator matrix
IMatrix              = ones(Num,Num);
for i = 1:1:Num
    Fi               = PopObj(i,:);
    Ir               = log(repmat(Fi,Num,1)./PopObj);
    MaxIr            = max(Ir,[],2);
    MinIr            = min(Ir,[],2);
    CVA              = MaxIr;
    DomInds          = find(MaxIr<=0);
    CVA(DomInds)     = MinIr(DomInds);
    IMatrix(:,i)     = CVA;
    IMatrix(i,i)     = Inf;
end

% calculate the size of individuals in the first nondominant level
IrFitness                = -min(IMatrix,[],2);
end
