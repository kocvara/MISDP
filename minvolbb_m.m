% Solve the multiple load minimum volume truss topology design problem with  binary
% variables using YALMIP's BNB solver
%
% Matlab coding by Michal Kocvara, University of Birmingham, 2010
% kocvara@maths.bham.ac.uk
%
% Input: structure "par" containing truss data; see the README file
%

m=par.m; n=par.n; n1=par.n1; BI=par.BI; xy=par.xy;
maska=par.maska; ijk=par.ijk;

% load case just given manually
ff=zeros(n1,1); ff1=ff; ff2=ff; ff3=ff; ff4=ff;
ff1(6)=-1;  ff2(12)=-1;
ff3(18)=-1; ff4(24)=-1;

compl = 10.0; par.cmp=compl; 

%t=sdpvar(m,1);
%t=binvar(m,1);
t=intvar(m,1);
assign(t,ones(m,1));

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

Constraints = [3 >= t >= 0];

Klarge1 = [compl -ff1'; -ff1 Kstiff];
Constraints = Constraints + [Klarge1 >= 0];

Klarge2 = [compl -ff2'; -ff2 Kstiff];
Constraints = Constraints + [Klarge2 >= 0];

Klarge3 = [compl -ff3'; -ff3 Kstiff];
Constraints = Constraints + [Klarge3 >= 0];

Klarge4 = [compl -ff4'; -ff4 Kstiff];
Constraints = Constraints + [Klarge4 >= 0];

Objective = t'*len;

Options=sdpsettings('solver','bnb','bnb.solver','mosek',...
    'bnb.maxiter',50000,...
    'usex0',1,'verbose',1,...
    'bnb.method','depthbest');  
%Options=sdpsettings('solver','mosek'); 

solvesdp(Constraints, Objective, Options);

t = double(t);

pic(par,t);


