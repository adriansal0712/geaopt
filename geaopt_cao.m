
% geaopt_cao: GEA optimization function (using Cao's method)
% [x_gea, f_gea, best_hist_gea] = geaopt_cao(pars,options)
% Enter parameters: 
% 1. pars: structure with compulsory fields fgname, nvar, extra elements as
% per the problem
% 2. options: structure with optional fields n (no of particles /
% solutions) and Num_iterations (no of iterations)
% Function should compute and generate the optimum value as the output [x_gea, f_gea, best_hist_gea]

function[x_gea, f_gea, best_hist_gea] = geaopt_cao(pars,options)        
 
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
    case {'F1','F2','F3','F4','F5','F6','F7','F8'}       
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

tic


% Uncomment/Comment the below line for Adrian's GEA Method
% [best_gea,fmin_gea,best_hist_gea] = gea_move_adr(Sol_gea,fmin_gea,fobj,best_gea,Fitness,lb,ub,n,Num_iterations,func_name,Q,x,y,z);

% uncomment/comment the below line for Cao's method
% Cao's method using frequency tuning for search, similar to Bat
    
[x_gea, f_gea, best_hist_gea] = gea_move_cao(pars,options,coordinates);

x_gea
f_gea

toc

end