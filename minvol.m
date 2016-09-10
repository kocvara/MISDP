% Solve the single load minimum volume truss topology design problem with
% continuous variables using YALMIP
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
compl = 1.0; par.cmp=compl; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t=sdpvar(m,1);
assign(t,(1.0/m).*ones(m,1));

len = zeros(m,1);
for i=1:m
   x1=xy(ijk(i,2)/2,1); y1=xy(ijk(i,2)/2,2);
   x2=xy(ijk(i,4)/2,1); y2=xy(ijk(i,4)/2,2);
   len(i,1)=sqrt((x1-x2)^2 + (y1-y2)^2);
end

%Ahelp=zeros(n,n);
%for i=1:m
%   Ahelp=Ahelp+len(i)*t(i)*BI(i,:)'*BI(i,:);
%end
%Astiff=Ahelp(maska,maska);
Astiff = BI'*diag(t.*len)*BI ;
Astiff=Astiff(maska,maska);

%F = set(1>t>0);

Alarge = [compl -ff'; -ff Astiff];
Constraints = [1>t>0, Alarge>0];
Objective = t'*len;

options=sdpsettings('solver','mosek',...
    'usex0',0,'verbose',1);

solvesdp(Constraints, Objective, options);

t = double(t);
%K = double(Astiff);
%ff'*(K\ff)

pic(par,t);


