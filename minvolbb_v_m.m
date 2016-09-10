% Solve the multiple mass minimum volume truss topology design problem with vibration constraints
% and  binary variables using YALMIP's BNB solver
%
% Matlab coding by Michal Kocvara, University of Birmingham, 2010
% kocvara@maths.bham.ac.uk
%
% Input: structure "par" containing truss data; see the README file
%

m=par.m; n=par.n; n1=par.n1;BI=par.BI;xy=par.xy;
maska=par.maska;ijk=par.ijk;

% PARAMETERS TO BE CHANGED MANUALLY (see also line 50 below)
compl = 0.1; par.cmp=compl; 
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

%ff=par.f; 
ff=zeros(n1,1);

mloc = [2 0 1 0; 0 2 0 1; 1 0 2 0; 0 1 0 2];
Massh = zeros(n,n);
for i=1:m 
    MA = sparse(n,n);
    mind = ijk(i,:);
    MA(mind,mind) = mloc*len(i);
    Massh = Massh + t(i)*MA;
end
Mass = Massh(maska,maska);

Mde1 = sparse(n1,n1);
Mde2 = sparse(n1,n1);

% PARAMETERS TO BE CHANGED MANUALLY
mm=300;
%Mde1(11,11) = mm; Mde1(12,12) = mm;Mde2(7,7) = mm; Mde2(8,8) = mm;
Mde1(17,17) = mm; Mde1(18,18) = mm;Mde2(23,23) = mm; Mde2(24,24) = mm;
%Mde1(31,31) = mm; Mde1(32,32) = mm;Mde2(39,39) = mm; Mde2(40,40) = mm;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eps=0.0;
F = [1>=t>=0];

Alarge = [compl -ff'; -ff Astiff];
F = [F; Alarge>0];

lambda_bar = 1.0e-3;

F = [F; Astiff - lambda_bar*(Mass+Mde1) >= 0];
F = [F; Astiff - lambda_bar*(Mass+Mde2) >= 0];
%F = [F; Astiff - lambda_bar*(Mass+Mde1+Mde2) >= 0];

options=sdpsettings('solver','bnb','bnb.solver','mosek',...
    'bnb.maxiter',50000,...
    'usex0',1,'verbose',1,...
    'bnb.method','depthbest');   

solvesdp(F,sum(t.*len),options);

t = double(t);
%K = double(Astiff);
%ff'*pinv(K)*ff

pic(par,t);


