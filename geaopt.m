
% geaopt: GEA optimization function
% [x_gea, f_gea, best_hist_gea] = geaopt(n,d,func_name,
% Num_iterations)
% Enter parameters: [n,d,func_name, Num_iterations]
% Function should compute and generate the optimum value as the output [best_gea, fmin_gea, best_hist_gea]

function[x_gea, f_gea, best_hist_gea] = geaopt(pars,options)        
 
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
    
[x_gea, f_gea, best_hist_gea] = gea_move_adr(pars,options,coordinates);

x_gea
f_gea

toc

end