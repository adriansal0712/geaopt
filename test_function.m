% Source: https://www.mathworks.com/matlabcentral/fileexchange/68981-bat-optimization-algorithm
% author: Abhishek Gupta
% modified by: Adrian Saldanha
% Function test_function is only for visualization of test functions. 
% For other optimization tasks, this function is invalid

function [range,dim,fobj,func_min] = test_function(F)

switch F
    case 'F1'               % De Jong's sphere function
        fobj = @F1;
%         lb = -100*ones(1,d);
%         ub = 100*ones(1,d);
        range = [-100 100 -100 100];
        dim=3;
%         y = x;
        func_min = 0;
        
    case 'F2'               % Schwefer 2.22 function
        fobj = @F2;
%         lb = -15*ones(1,d);
%         ub = 15*ones(1,d);
        range = [-15 15 -15 15];
        dim = 3;
%         y = x;
        func_min = 0;
        
    case 'F3'               % Griewangk's function
        fobj = @F3;
%         lb = -15*ones(1,d);
%         ub = 15*ones(1,d);
        range = [-15 15 -15 15];
        dim=3;
%         y = x;
        func_min = 0;
        
    case 'F4'               % Rosenbrock's function
        fobj = @F4;
%         lb = -15*ones(1,d);
%         ub = 15*ones(1,d);
        range = [-15 15 -15 15];
        dim = 3;
%         y = x;
        func_min = 0;
        
    case 'F5'               % Rastrigin's function
        fobj = @F5;
%         lb = -5*ones(1,d);
%         ub = 5*ones(1,d);
        range = [-5 5 -5 5];
        dim = 3;
%         y = x;
        func_min = 0;
        
    case 'F6'               % Michalewicz function
        fobj = @F6;
%         lb = 0*ones(1,d);
%         ub = 4*ones(1,d);
        range = [0 4 0 4];
        dim = 3;
%         y = x;
        func_min = -1.8013;
        
    case 'F7'
        
        fobj = @F7;
        range = [0 400 -40 40];

        dim = 3;
        
        func_min = -2.3291;
        
end
end 


function o = F1(x,y)          % De Jong's sphere function
o=x.^2 + y.^2;
end

function o = F2(x,y)          % Schwefer 2.22 function
o=abs(x)+abs(y)+prod(abs(x).*abs(y));
end

function o = F3(x,y)          % Griewangk's function
% dim = size(x,2);
dim = size(x,2);
o = 1+sum((x.^2)./4000) + sum((y.^2)./4000)-cos(x).*cos(y./sqrt(2));
end

function o = F4(x,y)          % Rosenbrock's function
dim=size(x,2);
o =((1-x).^2)+(100*((y-(x.^2)).^2));
end

function o = F5(x,y)          % Rastrigin's function
dim = size(x,2);
o = 10*2 + x.^2 + y.^2 - 10*cos(2*pi*x) - 10*cos(2*pi*y);
end

function o = F6(x, y)          % Michalewicz function
o = -sin(x).*(sin(x.^2./3.14159)).^20-sin(y).*(sin(2.*y.^2./3.14159)).^20;

end

function o = F7(x,y)
load example3_article.mat


C1=     [0 0 0 0 -1 0 1 0];
C2=     [0 0 0 0 1 0 0 0;
        0 0 0 0 0 1 0 0];

w = 19;                         % frequency of vibration to be suppressed

tau1 = 0;
tau2 = 0;

values = [w*j,-w*j];

for i = 1:size(x,1)
    for k = 1:size(y,1)
        
        K2 = [x(i,k),y(i,k)]';
        % tds = time-delay system with K2 as feedback
        tds = tds_create({A,B2*K2'*C2},[hA,tau2],{B},[0],{C},[0]);

        [g,tau1,R] = parameterizeDR(tds,values,B1,C1); % obtain 

        tds3 = tds_create({A,g*B1*C1,B2*K2'*C2},[hA tau1 tau2],{B},[0],{C},[0]);
        tds3.hD = [];
        options.minimal_real_part = -0.5;

        eigenvalues = compute_roots_DDAE(tds3,options);

        o(i,k) = max(real(eigenvalues.l0));

    end
end


end
    
