% Solve the single load minimum volume truss topology design problem with vibration constraints
% and  binary variables using YALMIP's BNB solver
%
% Matlab coding by Michal Kocvara, University of Birmingham, 2010
% kocvara@maths.bham.ac.uk
%
% Input: structure "par" containing truss data; see the README file
%
m=par.m; n=par.n; n1=par.n1; BI=par.BI; xy=par.xy;
maska=par.maska; ijk=par.ijk;

ff=par.f; 

% PARAMETERS TO BE CHANGED MANUALLY
compl = 2.0; par.cmp=compl; 
lambda_bar = 1.0e-2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t=binvar(m,1);
assign(t,ones(m,1));

len = zeros(m,1);
for i=1:m
   x1=xy(ijk(i,2)/2,1); y1=xy(ijk(i,2)/2,2);
   x2=xy(ijk(i,4)/2,1); y2=xy(ijk(i,4)/2,2);
   len(i,1)=sqrt((x1-x2)^2 + (y1-y2)^2);
end

Ahelp=zeros(n,n);
for i=1:m
   Ahelp=Ahelp+len(i)*t(i)*BI(i,:)'*BI(i,:);
end
Astiff=Ahelp(maska,maska);

mloc = [2 0 1 0; 0 2 0 1; 1 0 2 0; 0 1 0 2];
Massh = zeros(n,n);
for i=1:m  
    MA = sparse(n,n);
    mind = ijk(i,:);
    MA(mind,mind) = mloc*len(i);
    Massh = Massh + t(i)*MA;
end
Mass = Massh(maska,maska);

Mde = sparse(n1,n1);

F = [1 >= t >= 0];

Alarge = [compl -ff'; -ff Astiff];
F = [F; Alarge>=0]

F = [F; Astiff - lambda_bar*(Mass+Mde) >= 0];

options=sdpsettings('solver','bnb','bnb.solver','mosek',...
    'bnb.maxiter',50000,...
    'usex0',1,'verbose',1,...
    'bnb.method','depthbest'); 

solvesdp(F,sum(t.*len),options);

t = double(t);
%K = double(Astiff);
%ff'*pinv(K)*ff

pic(par,t);


