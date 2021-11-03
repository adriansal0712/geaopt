% psoopt: PSO optimization function
% [best_pso, fmin_pso, best_hist_pso] = psoopt(n,d,func_name,
% Num_iterations)
% Enter parameters: [n,d,func_name, Num_iterations]
% Function should compute and generate the optimum value as the output [best_pso, fmin_pso, best_hist_pso]

function [x_pso, f_pso, best_hist_pso] = psoopt(pars,options)

d = pars.nvar;
func_name = pars.fgname;

if ~isfield(options, 'Num_iterations')
    Num_iterations = 10*d;
    options.Num_iterations = Num_iterations;
else 
    Num_iterations = options.Num_iterations;
end

if ~isfield(options, 'n')
    n = 10*d;
    options.n = n;
else 
    n = options.n;
end

% Plotting function

switch func_name
    case {'F1','F2','F3','F4','F5','F6','F7'}
        
        figure
        [coordinates,range,func_min] = test_func_plot(func_name);    % Show function curve in 3-d

    otherwise
        coordinates.xplot = [];
        coordinates.yplot = [];
        coordinates.zplot = [];
        range = [];
        func_min = [];
        
end


% Initialization
tic
% Generating initial locations of n particles

% [Sol_pso, fitness_pso, fobj, lb,ub] = initpar(pars,options); % Call init_pso function to initialize solutions
% 
% fitness_best_sol = fitness_pso;         % Initialize historical best fitness of each particle
% Sol_best = Sol_pso;                     % Initialize historical best position
% 
% [f_pso, idx_pso] = min(fitness_pso);
% x_pso = Sol_pso(idx_pso,:);

% Plotting contour and visualization
   
% Start of Move function %

[x_pso, f_pso, best_hist_pso] = ...
        pso_move(pars,options,coordinates);

x_pso    
f_pso
    
toc
end
