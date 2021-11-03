% batopt: BAT optimization function
% [best_bat, fmin_bat, best_hist_bat] = batopt(n,d,func_name,
% Num_iterations)
% Enter parameters: [n,d,func_name, Num_iterations]
% Function should compute and generate the optimum value as the output [best_bat, fmin_bat, best_hist_bat]

function [x_bat, f_bat, best_hist_bat] = batopt(pars,options)

d = pars.nvar;
func_name = pars.fgname;

if ~isfield(options, 'Num_iterations')
    Num_iterations = 10*d;
else 
    Num_iterations = options.Num_iterations;
end

if ~isfield(options, 'n')
    n = 10*d;
else 
    n = options.n;
end

% Plotting function

switch func_name
    case {'F1','F2','F3','F4','F5','F6','F7'}       
        % Plot optimization curves, valid for predefined functions only
        
        figure
        [coordinates,range,func_min] = test_func_plot(func_name);    % Show function curve in 3-d

    otherwise
        
        % for other functions, generate blank coordinates
        
        coordinates.x = [];
        coordinates.y = [];
        coordinates.z = [];
        range = [];
        func_min = [];
        
end


% Generating initial locations of n particles
tic 

Q=zeros(n,1);                                   % frequency - for Cao's method only

[Sol_bat, fitness_bat, fobj, lb,ub] = initpar(pars,options);   % Call init_gea function to 
                                                % initialize solutions

% Find current best
[f_bat, I] = min(fitness_bat);
x_bat = Sol_bat(I,:);

% Plotting contour and visualization
   
% Start of Move function %
% BAT MOVE Function
[x_bat, f_bat, best_hist_bat] = bat_move(pars,options,coordinates);

x_bat    
f_bat

toc

end