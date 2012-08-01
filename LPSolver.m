function [x,V,Vbar,fval,exitflag,outputdata] = LPSolver(S0,t,T,threshold,eta,sim,type,Ppathselectedhour,Lpathselectedhour,Fselectedhour,YPathselectedhour)

B = 10000000;

% Choose what time points you want to sample, for the x and V values.
xselect = 1;
Vselect = 2;

% Preallocate vectors
a = zeros(T);
r = zeros(T);
y = zeros(T);
f = zeros(T);
P = zeros(T);
d = zeros(T);

% Process data from files
for k = 1:T
    a(k) = threshold;
    r(k) = Fselectedhour(k, sim);
    y(k) = YPathselectedhour(k, sim);
    f(k) = Fselectedhour(k, sim);
    d(k) = Lpathselectedhour(k, sim);
end
P(T) = Ppathselectedhour(type,sim);

% Compute YT, which is a sum of y(t)'s
YT = 0;
for j = 1:T-1
    YT = YT + y(j);
end

% Pre-allocate results vectors
x = zeros(T,1);
S = zeros(T,1);
V = zeros(T,1);
Vbar = zeros(T,1);

% Create vector of variables names
variables = cell(1,4*T);
N = length(variables);
for i = 1:T
    variables(1,i)     = {sprintf('x%u',i)};
    variables(1,T+i)   = {sprintf('S%u',i)};
    variables(1,2*T+i) = {sprintf('V%u',i)};
    variables(1,3*T+i) = {sprintf('Vbar%u',i)};
end

% For each variable get its name and position in the variables vector
for i = 1:N 
   eval([variables{i},' = ', num2str(i),';']); 
end

% Create the lower bound vector lb
 lb = zeros(size(variables));
 for i = 1:T
     % lb(1,i)     = 0;
     lb(1,T+i)     = -Inf;
     lb(1,2*T+i)   = -Inf;
     % lb(1,3*T+i) = 0;
 end
 
 % Create the upper bound vector as a vector of Inf
 ub = Inf(size(variables));
 % Add the x's upper bounds that are needed to avoid infinite solutions.
 for i = 1:T
     ub(1,i) = B; 
 end
 
 % Create the linear inequalities matrix A and constants vector b such that A*Vars<=b.
 % In the LP formulation there are only T of these inequality constraints.
 A = zeros(T,N);
 b = zeros(T,1);
 for i = 1:T-1
     eval(sprintf('A(%u,V%u)=1; A(%u,Vbar%u)=-1;',i,i,i,i));
     b(i,1) = a(i);
 end
 eval(sprintf('A(T,Vbar%u)=-1; A(T,S%u)=-P(T);',T,T));
 b(T,1) = a(T)-P(T)*(d(T)+YT);
 
 % Create the linear equalities matrix Aeq and constants vector beq such that Aeq*Vars=beq.
 % In the LP formulation there are 2*T+1 of these equality constraints.
 % We will code them in two groups plus one special single constraint.
 Aeq = zeros(2*T+1,N);
 beq = zeros(2*T+1,1);
 
 % First group of equality constraints: fx constraints. There are T of these.
 for i = 1:T-1
     eval(sprintf(...
         'Aeq(%u,x%u)=f(i); Aeq(%u,V%u)=-1; Aeq(%u,V%u)=1; Aeq(%u,Vbar%u)=eta;',...
         i,i,i,i,i,i+1,i,i+1));
     beq(i,1) = r(i)*y(i);
 end
 eval(sprintf('Aeq(T,V%u)=1; Aeq(T,S%u)=P(T); Aeq(T,Vbar%u)=-eta;',T,T,T));
 beq(T,1) = P(T)*(d(T)+YT);
 
 % Second group of equality constraints: Sk constraints. There are T of these.
 for i = 1:T
     eval(sprintf('Aeq(%u,S%u)=1; Aeq(%u,1:%u)=-1;',T+i,i,T+i,i-1));
     beq(T+i,1) = S0;
 end
 
 % Special xT=0 constraint.
 eval(sprintf('Aeq(%u,x%u)=1;',2*T+1,T)); 
 % beq(1,2*T+1)=0;
 
 % Define the objective function vector g such that we minimize g*Vars.
 g = zeros(size(variables));
 g(1,V1) = 1;
 
 % Solve the LP
 %[Opt_Sol fval exitflag] = linprog(g,A,b,Aeq,beq,lb,ub);
 [Opt_Sol fval exitflag] = cplexlp(g,A,b,Aeq,beq,lb,ub);

% Extract results from Opt_Sol and put into results variables
for i = 1:T
    x(i,1) = Opt_Sol(i);
    S(i,1) = Opt_Sol(T+i);  
    V(i,1) = Opt_Sol(2*T+i);
    Vbar(i,1) = Opt_Sol(3*T+i);
end

% Fetch x(xselect) and V(Vselect), xselect and Vselect to be determined by the user
outputdata = [x(xselect);V(Vselect)];

end


