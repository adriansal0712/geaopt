function [fitness] = relaysim_gea(Sol, i)

        simOut = sim('relay_model_gea_2017a.slx', 'ReturnWorkspaceOutputs', 'on');
        fitness = simOut.err_integral.Data(end);    

end