function o = relaysim(x,pars)
% Function relaysim: Function to simulate Simulink model of relay feedback
%   Detailed explanation goes here

addpath("relay_id");

% PID Parameters - initialization
assignin('caller','rp',1);
assignin('caller','ri',0);
assignin('caller','rd',0);

% Step input to system - initialization
assignin('caller','step_in',0);
         

handle = load_system('relay_model_gea_2017a');
set_param('relay_model_gea_2017a/sw1', 'sw', '0');
if isfield(pars,'process')
    
    global s;
    syms s;
    process = pars.process;
    [G,relay,Tds,num,den,Tsim,tol] = get_function(process,s);
    
end

K = x(1);
a2 = x(2);
a1 = x(3);
Td = x(4);

simOut = sim('relay_model_gea_2017a.slx','ReturnWorkspaceOutputs','on');
o = simOut.err_integral.Data(end);

end

