function cns(param_name)
    close all
    rand('state',sum(100*clock));
    param = param_initial(param_name);
    
    if strcmp(param_name,'ss00')==1
        neuron(param);
    elseif strcmp(param_name,'network')==1
		neuron_network(param);
    elseif strcmp(param_name,'pattern_positive')==1
    
    elseif strcmp(param_name,'pattern_negative')==1
        
    else
        fprintf('Wrong Argument\n');
        return;
    end
end
