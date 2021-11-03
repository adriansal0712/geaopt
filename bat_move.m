% References: ----------------------------------------------------------- %
% 1) Yang X.S. (2010). A New Metaheuristic Bat-Inspired Algorithm, In:    %
%   Nature-Inspired Cooperative Strategies for Optimization (NICSO 2010), %
%   Studies in Computational Intelligence, vol 284. Springer, Berlin.     %
%   https://doi.org/10.1007/978-3-642-12538-6_6                           %
% 2) Yang X.S. (2014). Nature-Inspired Optimization Algorithms, Elsevier. %
% ----------------------------------------------------------------------- %
% move all particles toward global minimum
function [x_bat, f_bat, best_hist] = bat_move(pars,options,coordinates)
                                           
                                         
% Obtain coordinates

if nargin<3
    coordinates.x = [];
    coordinates.y = [];
    coordinates.z = [];
end

xplot = coordinates.x;
yplot = coordinates.y;
zplot = coordinates.z;

figure

Num_iterations = options.Num_iterations;    % number of iterations
Tmax = Num_iterations;
n = options.n;                              % no of solutions
func_name = pars.fgname;                    % func_name
d = pars.nvar;                              % nvar of pars 

% Parameterization
ro = 1;                         % pulse-emission rate
A = 1;
Qmax = 2;               % Max frequency
Qmin = 0;               % Min frequency


[Sol, Fitness, fobj, lb,ub] = initpar(pars,options);   % Call initpar function to 
                                                % initialize solutions

alpha = 0.97;              % constant for loudness update = 0.97
gamma = 0.1;              % constant for emission rate update
sigma = 0.1*(ub - lb);                % local search scope
theta = 0.7;

% Initialization
Q = zeros(n,1);                   % initialize frequencies
v = zeros(n,d);                   % initialize velocities

% Find current best
[f_bat, I] = min(Fitness);
x_bat = Sol(I,:);


for iter = 1: Num_iterations,
    
    switch func_name
        case 'F1' 
            
            hold off;
            contour(xplot,yplot,zplot, 15); 
            hold on;
            title('De Jong sphere function using BAT', 'FontSize', 16); 
            plot(Sol(:,1), Sol(:,2), 'b.', x_bat(1), x_bat(2),'r*');
            axis([lb(1) ub(1) lb(2) ub(2)]);


        case 'F2'
            hold off;
            contour(xplot,yplot,zplot, 15); 
            hold on;
            title('Schwefel 2.2 function using BAT', 'FontSize', 16); 
            plot(Sol(:,1), Sol(:,2), 'b.', x_bat(1), x_bat(2),'r*');
            axis([lb(1) ub(1) lb(2) ub(2)]);

        case 'F3' 
            hold off;
            contour(xplot,yplot,zplot, 15); 
            hold on;
            title('Griewangk function using BAT', 'FontSize', 16); 
            plot(Sol(:,1), Sol(:,2), 'b.', x_bat(1), x_bat(2),'r*');
            axis([lb(1) ub(1) lb(2) ub(2)]);

        case 'F4'
            hold off;
            contour(xplot,yplot,zplot, 15); 
            hold on;
            title('Rosenbrock function using BAT', 'FontSize', 16); 
            plot(Sol(:,1), Sol(:,2), 'b.', x_bat(1), x_bat(2),'r*');
            axis([lb(1) ub(1) lb(2) ub(2)]);

        case 'F5' 
            hold off;
            contour(xplot,yplot,zplot, 15); 
            hold on;
            title('Rastrigin function using BAT', 'FontSize', 16); 
            plot(Sol(:,1), Sol(:,2), 'b.', x_bat(1), x_bat(2),'r*');
            axis([lb(1) ub(1) lb(2) ub(2)]);

        case 'F6'
            hold off;
            contour(xplot,yplot,zplot, 15); 
            hold on;
            title('Michalewicz function using BAT', 'FontSize', 16); 
            plot(Sol(:,1), Sol(:,2), 'b.', x_bat(1), x_bat(2),'r*');
            axis([lb(1) ub(1) lb(2) ub(2)]);
            
        case 'F7'
            hold off;
            contour(xplot,yplot,zplot, 50); 
            hold on;
            title('Michalewicz function using BAT', 'FontSize', 16); 
            plot(Sol(:,1), Sol(:,2), 'b.', x_bat(1), x_bat(2),'r*');
            axis([lb(1) ub(1) lb(2) ub(2)]);
            
        case 'F8'
            hold off;
            contour(xplot,yplot,zplot, 50); 
            hold on;
            title('Non-collocated control using GEA', 'FontSize', 16); 
            plot(Sol(:,1), Sol(:,2), 'b.', x_bat(1), x_bat(2),'r*');
            
        otherwise
            scatter(iter,f_bat,'b.');
            hold on
            title('Cost Function');
            xlim([0 Tmax]);
            xlabel("Number of iterations");
            ylabel("Cost Function");
            
    end
    
    S = Sol;
    
    A = alpha*A;                                % updating the loudness
    rho = ro*(1-exp(-gamma*iter));              % updating the pulse emission rate

    Q = Qmin + (Qmax-Qmin)*rand(n,1);           % update frequencies randomly
    v = v.*theta + (Sol - x_bat).*Q;     % update velocities
    
    % Simple bounds
    Flag4up=S>ub;
    Flag4low=S<lb;
    S=S.*(~(Flag4up+Flag4low))+ub.*Flag4up+lb.*Flag4low;
        

    for i = 1:n

        if rand < rho,
            S(i,:) = x_bat + sigma.*randn(1,d)*A;    % local search around best
        end
        
     % Simple bounds                                 % check if positions exceed design space
        Flag4up = S(i,:) > ub;                
        Flag4low = S(i,:) <lb;
        S(i,:) = S(i,:).*(~(Flag4up+Flag4low)) + ub.*Flag4up + lb.*Flag4low;

    % 3 - Generate new solution by flying randomly
    % 
    %         epsi = unifrnd(-1,1, [1 2]);
    %         S_bat(i,[1 2]) = Sol(i,[1 2]) + epsi*A;  % *mean(A);

        fitnessnew = fobj(S(i,:),pars); % calculate fitness function

        if fitnessnew < Fitness(i) & rand>A,
            Fitness(i) = fitnessnew;
            Sol(i,:) = S(i,:);
        end

        % if minimum value is lower than the global minimum, select individual
        if Fitness(i) <= f_bat,
            x_bat = Sol(i,:);
            f_bat = Fitness(i);
        end
            
    end


    drawnow;

    best_hist(iter,:) = [iter, f_bat];

end