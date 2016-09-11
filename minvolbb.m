% Solve the single load minimum volume truss topology design problem with  binary
% variables using YALMIP's BNB solver
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
compl = 2*6.1; compl = 1; par.cmp=compl; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t=intvar(m,1);
assign(t,ones(m,1));;

len = zeros(m,1);
for i=1:m
   x1=xy(ijk(i,2)/2,1); y1=xy(ijk(i,2)/2,2);
   x2=xy(ijk(i,4)/2,1); y2=xy(ijk(i,4)/2,2);
   len(i,1)=sqrt((x1-x2)^2 + (y1-y2)^2);
end

Khelp=zeros(n,n);
for i=1:m
   Khelp=Khelp+len(i)*t(i)*BI(i,:)'*BI(i,:);
end
Kstiff=Khelp(maska,maska);

F = [1 >= t >= 0];

Klarge = [compl -ff'; -ff Kstiff];
F = [F; Klarge>=0];

options=sdpsettings('solver','bnb','bnb.solver','mosek',...
    'bnb.maxiter',50000,...
    'usex0',1,'verbose',1,...
    'bnb.method','depthbest');

solvesdp(F,sum(t.*len),options);

t = double(t);
%K = double(Kstiff);
%ff'*pinv(K)*ff;

pic(par,t);



