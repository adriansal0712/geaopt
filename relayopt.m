

function [best, fmin, best_hist] = relayopt(n,d,func_name, Num_iterations)

% Obtain function parameters depending on function name

syms s;
[Gs, relay, Tds, num, den, Tsim, tol] = get_function(func_name, s);

% Sending values to base workspace
assignin('base', 'num', num);
assignin('base', 'den', den);
assignin('base', 'Tds', Tds);
assignin('base', 'Gs', G);
assignin('base', 'relay', relay);
assignin('base', 'Tsim', Tsim);
assignin('base', 'tol', tol);



end
