% Author: Adrian Saldanha
% Objective function: For obtaining the default parameters based on the
% function name.
% Accepts inputs: (x,d), where x = func_name, d = no. of dimensions
% Assigns outputs: [lb, ub, fobj, func_min]

function [lb,ub,fobj,func_min] = fobjective(pars)

d = pars.nvar;
func_name = pars.fgname;

switch func_name
    case 'F1'                   % De-Jong's Sphere Function
        fobj = @F1;
        lb = -100*ones(1,d);
        ub = 100*ones(1,d);
        func_min = 0;
        
    case 'F2'               % Schwefer 2.22 function
        fobj = @F2;
        lb = -15*ones(1,d);
        ub = 15*ones(1,d);
        func_min = 0;
        
    case 'F3'               % Griewangk's function
        fobj = @F3;
        lb = -15*ones(1,d);
        ub = 15*ones(1,d);
        func_min = 0;
        
    case 'F4'               % Rosenbrock's function
        fobj = @F4;
        lb = -15*ones(1,d);
        ub = 15*ones(1,d);
        func_min = 0;
        
    case 'F5'               % Rastrigin's function
        fobj = @F5;
        lb = -5*ones(1,d);
        ub = 5*ones(1,d);
        func_min = 0;
        
    case 'F6'               % Michalewicz function
        fobj = @F6;
        lb = 0*ones(1,d);
        ub = 4*ones(1,d);
        func_min = -1.8013;
        
    case 'F7'
        
        fobj = @F7;
        range = [0 400 -40 40];
        lb = [range(1),range(3)];
        ub = [range(2),range(4)];

        func_min = -2.3291;
        
    case 'F8'
        
        fobj = @F8;
        range = [0 100 -100 0];
        lb = [range(1),range(3)];
        ub = [range(2),range(4)];

        func_min = -0.3943;
        
    otherwise
        fobj = @FX;
        lb = -1000*ones(1,d);
        ub = 1000*ones(1,d);
        
        func_min = [];
        
end
end

function o = F1(x,pars)          % De Jong's sphere function
    o=sum(x.^2);
end

function o = F2(x,pars)          % Schwefer 2.22 function
o=sum(abs(x))+prod(abs(x));
end

function o = F3(x,pars)          % Griewangk's function

s = 0;
p = 1;
d = size(x,2);

for i=1:d,
    s = s+sum(x(i)^2/4000);
    p = p*cos(x(i)/sqrt(i));
    
end
o = 1+s-p;

end

function o = F4(x,pars)          % Rosenbrock's function
d=size(x,2);

s = 0; p = 1;

for i = 1:d-1,
    s = s+sum(100*(x(i+1) - x(i)^2)^2 + (x(i) - 1)^2);
end

o = s;
end

function o = F5(x,pars)          % Rastrigin's function
d = size(x,2);

s = 0;
for i = 1:d,
    s = s + x(i)^2 - 10*cos(2*pi*x(i));
end

o = 10*d + s;

end

function o = F6(x,pars)          % Michalewicz function

d = size(x,2);
s = 0;

for i = 1:d,
    s = s-sin(x(i))*(sin((i*(x(i))^2)/pi))^20;
end
    
o = s;

end

function o = F7(x,pars)
load example3_article.mat


C1=     [0 0 0 0 -1 0 1 0];
C2=     [0 0 0 0 1 0 0 0;
        0 0 0 0 0 1 0 0];

w = 19;                         % frequency of vibration to be suppressed

tau1 = 0;
tau2 = 0;

values = [w*j,-w*j];


K2 = x';
% tds = time-delay system with K2 as feedback
tds = tds_create({A,B2*K2'*C2},[hA,tau2],{B},[0],{C},[0]);

[g,tau1,R] = parameterizeDR(tds,values,B1,C1); % obtain 

tds3 = tds_create({A,g*B1*C1,B2*K2'*C2},[hA tau1 tau2],{B},[0],{C},[0]);
tds3.hD = [];
options.minimal_real_part = -0.5;

eigenvalues = compute_roots_DDAE(tds3,options);

o = max(real(eigenvalues.l0));


end

function o = F8(x,pars)

    options = tdsrootsoptions;                  % Options for computing roots
    options.minimal_real_part = -0.2;   
    options.fixN=0;
    options.minimal_real_part = -0.5;

    o = compute_gradient(x,pars);


end

function o = FX(x,pars)

    x = x';
    func = str2func(pars.fgname);
    o = func(x,pars);

end
