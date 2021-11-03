%% Init PAR function - for initializing the parameters of GEA/BAT/PSO algorithm
% Standard procedure for initializing, applicable for all three algorithms
% - GEA, BAT, PSO

function [Sol, Fitness, fobj, lb,ub] = initpar(pars,options);

n = options.n;
d = pars.nvar;
func = pars.fgname;

[lb,ub,fobj,func_min] = fobjective(pars);

Sol = rand(n,d).*(ub-lb) + lb;

Fitness = zeros(n,1);

for i = 1:n,
    Fitness(i,1) = fobj(Sol(i,:),pars);     %  calculating fitness
end


end
