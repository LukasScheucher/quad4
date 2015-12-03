%% SETUP 1
%       
%
%  ^y
%  |     Master
%  |
% (1)-----(2)-----(3)         
%  |       |       |          %%  ^y
%  |  [1]  |  [2]  |          %%  |     Slave
%  |       |       |          %%  |
% (4)-----(5)-----(6)         %% (o)----(o)-----(o)
%  |       |       |          %%  |       |      |
%  |  [3]  |  [4]  |          %%  |  [X]  | [X]  |
%  |       |       |          %%  |       |      |
% (7)-----(8)-----(9)         %% (o)----(o)-----(o)
%  |       |       |          %%  |       |      |
%  |  [5]  |  [6]  |          %%  |  [X]  | [X]  |
%  |       |       |          %%  |       |      |
% (10)-----(11)-----(12)-->x  %% (o)----(o)-----(o)----->x
%

sizefac=2; %determines, how many times smaller the slave elements are

nodes=[...
0 2 4 0  2 4 0 2  4 0 2 4
6 6 6 4  4 4 2 2  2 0 0 0]';

elements=[...
    4  5  2  1
    5  6  3  2
    7  8  5  4
    8  9  6  5
    10 11 8  7
    11 12 9  8];


numslaveele=sizefac*sizefac*4;



% o----------->i
% |
% |
% |
% \/
% j


%% creating nodes
for j=1:2*sizefac+1
  for i=1:2*sizefac+1
     nodes=[...
            nodes
            (i-1)*2/sizefac, 4-(j-1)*2/sizefac ];
     
  end
end



%% creating elements
temp=sizefac*2+1;
for j=1:2*sizefac
  for i=1:2*sizefac
    elements=[...
      elements
      [j*temp+1*i j*temp+i+1 (j-1)*temp+i+1 (j-1)*temp+i]+12 ];
  end
end


numnodes=length(nodes);
numele=length(elements);

% params.bodyele  =[0 0 2 2 2 2 1 1 1 1];
% params.ifaceele =[0 0 2 2 2 2 1 1 1 1];
% params.bodynode =[2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1];
% params.ifacenode=[0 0 0 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1];

params.bodyele  =[ones(1,6)*2 ones(1,sizefac*sizefac*4)];
params.ifaceele =[zeros(1,2) 2*ones(1,4) ones(1,sizefac*sizefac*4)];

params.bodynode = [2*ones(1,12) ones(1,numnodes-12)];
params.ifacenode= [zeros(1,3) 2*ones(1,9) ones(1,numnodes-12)];

params.dofshift=(2*sizefac+1)^2*2;
params.stresstype='planestress',
params.E=1000,
params.v=0.3,

params.alpha=@(y) 0.5;%y./4; %auf master bezogen



dirichdof=[...
           dofs(find(nodes(:,2)==6)) [dofs(find(nodes(:,2)==6))*0]
           dofs(numnodes-sizefac*2:numnodes) 0.2*ones(length(dofs(numnodes-sizefac*2:numnodes)),1)];

       
%% modification

% No modification
% params.modnodes=[]
% params.modele =[]
% params.ifacenode(params.modnodes)=3
% params.dofshift=params.dofshift



params.modnodes=numnodes-sizefac*2:numnodes
params.modele =numele-sizefac*2+1:numele
params.ifacenode(params.modnodes)=3
params.ifaceele(params.modele)=3;

         
% disp=zeros(length(nodes),2);
% f=figure();
% plotdisp(f,params,elements,nodes,disp)
       
 