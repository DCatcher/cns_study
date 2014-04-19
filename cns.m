function cns(param_name)
    close all
    %rand('state',sum(100*clock));
    rand('state',0);
    param = param_initial(param_name);
    
    if strcmp(param_name,'ss00')==1
        neuron(param);
    elseif strcmp(param_name,'ss01')==1
%		neuron_network(param);
		neuron_ss01(param);
    elseif strcmp(param_name,'pattern_positive')==1
		neuron_pattern_positive(param);
	elseif strcmp(param_name,'sparse_coding')==1
		neuron_sparse_coding(param);
    elseif strcmp(param_name,'pattern_negative')==1
		fprintf('I will finish it!');
		return;
    else
        fprintf('Wrong Argument\n');
        return;
    end
end
