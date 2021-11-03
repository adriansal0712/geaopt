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
        
    case 'F8'
        
        fobj = @F8;
        range = [0 500 -500 100];

        dim = 3;
        
        func_min = -0.3943;
        
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


function o = F8(x,y)

% Haik's Setup
    % Mass
    ma = 0.5; m0 = 1; m1 = 0.5; m2 = 1.15;

    % Stiffness
    ka = 400; k0 = 410; k1 = 1450; k2 = 380; k3 = 405; kp = 1500;

    % Damping
    ca = 1.9; c0 = 2.1; c1 = 4.9; c2 = 2.2; c3 = 2; cp = 5;

    % State-matrix
    A = [0  1   0   0   0   0   0   0;
        -(k0+k3+kp)/m0  -(c0+c3+cp)/m0  k0/m0   c0/m0   k3/m0   c3/m0   0   0;
        0   0   0   1   0   0   0   0;
        k0/m1   c0/m1   -(k0+k1)/m1     -(c0+c1)/m1     k1/m1   c1/m1   0   0;
        0   0   0   0   0   1   0   0;
        k3/m2   c3/m2   k1/m2   c1/m2   -(k1+k2+k3+ka)/m2   -(c1+c2+c3+ca)/m2   ka/m2   ca/m2;
        0   0   0   0   0   0   0   1;
        0   0   0   0   ka/ma   ca/ma   -ka/ma  -ca/ma]


    % System input matrix
    B = [0  1/m0    0   0   0   0   0   0]';
    % System output matrix
    C = [0  0   1   0   0   0   0   0];

    % Feed-through matrix
    D = 0;

    % Input Matrix (Controllability)
    B1 = [
            0   0   0   0   0  -1/m2 0   1/ma;
            ]';
    C1 = [  
            0   0   0   0   -1  0   1   0;
            0   0   0   0   1   0   0   0;
            0   0   0   0   0   1   0   0;
            1   0   0   1   0   0   0   0
            ];


    n = size(A,1);  % order of the system
    p = size(C1,1); % degrees of freedom of controller = size(K,2)
    q = size(B1,2); % no of control inputs
    hA = [0];

    tds = tds_create({A},hA,{B},[0],{C},[0]);
    options = tdsrootsoptions;                  % Options for computing roots
    options.minimal_real_part = -0.2;   
    options.fixN=0;
    tds_check_valid(tds);

    w = 19;         % frequency of vibration ot be supressed
    values = sort([w*j,-w*j]);
    m = size(values,2);     % number of zeros to be placed

    ipdelay = 0;

    pars = parInput(tds,B1,C1,ipdelay,w);
              
    
    for i = 1:size(x,1)
        for k = 1:size(y,1)

            % Reshaping K2
            K2 = zeros(p,q);                    % matrix of dim(pxq) = K2'
            K2(m+1:end) = [x(i,k),y(i,k)]';     % elements 1:m are zeros
                                                % elements m+1:end take the values from opt
            K2 = K2';
    
            tds2 = tds_create({tds.A{:},B1*K2*C1},[tds.hA,ipdelay],...
                    {tds.B1{:}},[tds.hB1],{tds.C1{:}},[tds.hC1]);
    
            [K0 G] = assign_zeros(tds2,values,B1,C1,ipdelay);
            K1 = [K0' zeros(1,(p-m))];      % [g^T 0]
            I = eye(q);                     % Identity
            e1 = I(:,1);                    % Unit vector to extract the 1st column of the B1

            K = e1*K1 + K2;      % Controller gain matrix

            tds3 = tds_create({tds.A{:},B1*K*C1},[tds.hA,ipdelay],{tds.B1{:}},[0],{tds.C1{:}},[0]);
            tds3.hD = [];
            options.minimal_real_part = -0.5;

            eigenvalues = compute_roots_DDAE(tds3,options);
            eigenvalues.l0;
            o(i,k) = max(real(eigenvalues.l0));
            
        end
    end

end
    
