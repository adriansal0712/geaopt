%% PSO Move function - for moving the particles to the new location
% move all particles toward 'best'

function [x_pso, f_pso, best_hist_pso] = ...
        pso_move(pars,options,coordinates)

    
% Obtain coordinates
if nargin<3
    coordinates.x = [];
    coordinates.y = [];
    coordinates.z = [];
end

xplot = coordinates.x;
yplot = coordinates.y;
zplot = coordinates.z;


% Parameterization - PSO
beta = 2;                         % medium speed of convergence
gamma = 0.97;                       % gamme: 0.9 -> 0.97
alpha0 = 2;
theta = 0.7;                        % inertia weight

figure

Num_iterations = options.Num_iterations;        % number of iterations
Tmax = Num_iterations;
n = options.n;                                  % no of solutions
func_name = pars.fgname; 
d = pars.nvar;                


[Sol, Fitness, fobj, lb,ub] = initpar(pars,options);   % Call init_PSO function to 

% Initialize best position and fitness for each particle
fitness_best_sol = Fitness;        
Sol_best = Sol;                     

% Find current best
[f_pso, I] = min(Fitness);
x_pso = Sol(I,:);

v = zeros(n,d);                     % 2 = 2x dimensions


for iter = 1: Num_iterations,

    %   % show contour of function
    switch func_name
        case 'F1' 
            
            hold off;
            contour(xplot,yplot,zplot, 15); 
            hold on;
            title('De Jong sphere function using PSO', 'FontSize', 16); 
            plot(Sol(:,1), Sol(:,2), 'b.', x_pso(1), x_pso(2),'r*');
            axis([lb(1) ub(1) lb(2) ub(2)]);


        case 'F2'
            hold off;
            contour(xplot,yplot,zplot, 15); 
            hold on;
            title('Schwefel 2.2 function using PSO', 'FontSize', 16); 
            plot(Sol(:,1), Sol(:,2), 'b.', x_pso(1), x_pso(2),'r*');
            axis([lb(1) ub(1) lb(2) ub(2)]);

        case 'F3' 
            hold off;
            contour(xplot,yplot,zplot, 15); 
            hold on;
            title('Griewangk function using PSO', 'FontSize', 16); 
            plot(Sol(:,1), Sol(:,2), 'b.', x_pso(1), x_pso(2),'r*');
            axis([lb(1) ub(1) lb(2) ub(2)]);

        case 'F4'
            hold off;
            contour(xplot,yplot,zplot, 15); 
            hold on;
            title('Rosenbrock function using PSO', 'FontSize', 16); 
            plot(Sol(:,1), Sol(:,2), 'b.', x_pso(1), x_pso(2),'r*');
            axis([lb(1) ub(1) lb(2) ub(2)]);

        case 'F5' 
            hold off;
            contour(xplot,yplot,zplot, 15); 
            hold on;
            title('Rastrigin function using PSO', 'FontSize', 16); 
            plot(Sol(:,1), Sol(:,2), 'b.', x_pso(1), x_pso(2),'r*');
            axis([lb(1) ub(1) lb(2) ub(2)]);

        case 'F6'
            hold off;
            contour(xplot,yplot,zplot, 15); 
            hold on;
            title('Michalewicz function using PSO', 'FontSize', 16); 
            plot(Sol(:,1), Sol(:,2), 'b.', x_pso(1), x_pso(2),'r*');
            axis([lb(1) ub(1) lb(2) ub(2)]);
            
        case 'F7'
            hold off;
            contour(xplot,yplot,zplot, 50); 
            hold on;
            title('Optimizing spectral abscissa using PSO', 'FontSize', 16); 
            plot(Sol(:,1), Sol(:,2), 'b.', x_pso(1), x_pso(2),'r*');
            axis([lb(1) ub(1) lb(2) ub(2)]);
            
        case 'F8'
            hold off;
            contour(xplot,yplot,zplot, 50); 
            hold on;
            title('Non-collocated control using GEA', 'FontSize', 16); 
            plot(Sol(:,1), Sol(:,2), 'b.', x_pso(1), x_pso(2),'r*');
            axis([lb(1) ub(1) lb(2) ub(2)]);
            
        otherwise
            scatter(iter,f_pso,'b.');
            hold on
            title('Cost Function');
            xlim([0 Tmax]);
            xlabel("Number of iterations");
            ylabel("Cost Function");
            
    end




    % 1 - Update velocities

    v(:,[1:d]) = theta*v(:,[1:d]) + alpha0*rand(n,d).*(x_pso - Sol(:,[1:d])) ...
        + beta*rand(n,d).*(Sol_best(:,[1:d]) - Sol(:,[1:d]));


    % 2 - Update positions

    S(:,[1:d]) = Sol(:,[1:d]) + v(:,[1:d]);


    % Applying simple bounds
    Flag4up=S(:,[1:d])>ub;
    Flag4low=S(:,[1:d])<lb;
    Sol(:,[1:d])=S(:,[1:d]).*(~(Flag4up+Flag4low))+ub.*Flag4up+lb.*Flag4low;


    % 3 - Evaluate objective function at new location

    for ii = 1:n
        fitness_pso(ii) = fobj(Sol(ii,[1:d]),pars);

        % 4 - Find current best for each particle
        
        if fitness_pso(ii) < fitness_best_sol(ii)
            fitness_best_sol(ii) = fitness_pso(ii);
            Sol_best(ii, :) = Sol(ii, :);
        end

        % 5 - Update global minimum

        if fitness_pso(ii) < f_pso,
            f_pso = fitness_pso(ii);
            x_pso = Sol(ii, :);
        end



    end
    
    drawnow;
    best_hist_pso(iter,:) = [iter, f_pso];
    


end