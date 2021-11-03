%% GEA Algorithm - MAIN FUNCTION
% Program revised with the following new aspects:
% 1: Functions work for multiple dimensions (not restricted to 2
% 2: Main function is clean - All operations take place nested in the
% sub-fuctions geaopt, batopt and psoopt
%% Initializations

clear,close all

%% Initialize general parameters

d = 10;                        % number of dimensions
n = 10*d;                   % number of solutions, for faster convergence, set n = 50
Num_iterations = 20*d;      % number of iterations


%% Plot function curve - for display only, in 3-d
% F1 = De Jong's sphere function
% F2 = Schwefel 2.2 function
% F3 = Griewangk's function
% F4 = Rosenbrock's function
% F5 = Rastrigin's function
% F6 = Michalewicz function


% fgname = 'P1';           % Enter function name F1...F6


% For system identification
% P1 = 1/((s+1)*(0.1*s+1)*(0.01*s+1)*(0.001*s+1))
% P2 = 1/(s+1)^4
% P3 = e^-s/(0.05*s+1)^2
% P4 = 1/0.5*s^2 + s + 1)^2
% P5 = (e^-s)/((s+1)*(0.3s+1)^2))
% P6 = (e^-4s)/((s+1)*(0.3*s+1)^2)

d = 4;                      % dimension of the search variable is 4
fgname = "relaysim";        % file name of optimization cost function
process = 'P1';             % For predefined processes
pars.process = process;

%% Inputting parameters to variables

pars.fgname = fgname;
pars.nvar = d;

options.n = 10*d;
options.Num_iterations = 20*d;


%% PART 1 - GUIDING EVOLUTIONARY ALGORITHM
% a) Using Modified GEA method - default
% b) Using Cao's Method - Go to function geaopt, comment out gea_move_adr
% and uncomment gea_move_cao


% Using the modified method
% best_gea = zeros(10,1);
% for i = 1:10,
    [x_gea,f_gea,best_hist_gea] = geaopt(pars,options);
%     best_gea(i) = f_gea;
% end
% avg_gea = mean(best_gea);
% 
% % Using Cao's Method:
% best_gea_cao = zeros(10,1);
% for i = 1:10,
    [x_gea_cao,f_gea_cao,best_hist_gea_cao] = geaopt_cao(pars,options);
%     best_gea_cao(i) = f_gea_cao;
% end
% avg_gea_cao = mean(best_gea_cao);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF GEA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PART 2 - Using Accelerated PSO Algorithm for comparison

% Setting the parameters: alpha, beta
% Random amplitude of roaming particles alpha = [0 1]
% alpha = gamma^t = 0.7^t
% Speed of convergence (0->1) = (slow -> fast)
   
% best_pso = zeros(10,1);
% for i = 1:10,
    [x_pso, f_pso, best_hist_pso] = psoopt(pars,options);
%     best_pso(i) = f_pso;
% end
% avg_best_pso = mean(best_pso);


%%%%%%%%%%%%%%%%%%%%%%%%%END OF PSO ALGORITHM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PART 3 - Using Bat Algorithm

% Start iterations

% best_bat = zeros(10,1);
% for i = 1:10,
    [x_bat, f_bat, best_hist_bat] = batopt(pars,options);
%     best_bat(i) = f_bat;
% end                                    
% avg_best_bat = mean(best_bat);

%%%%%%%%%%%%%%%%%%%%%%%%END OF BAT ALGORITHM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot graph of minima vs iteration number

figure
plot_iter_curve(best_hist_gea, best_hist_pso, best_hist_bat, Num_iterations, fgname);

text = "results/"+fgname+".mat";

save(text,"best_gea","avg_gea","best_gea_cao","avg_gea_cao"...
    ,"best_pso","avg_best_pso","best_bat","avg_best_bat");


