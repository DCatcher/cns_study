function cns(param_name)
    close all
    %rand('state',sum(100*clock));
    rand('state',0);
	randn('state',0);
    param = param_initial(param_name);
    
    if strcmp(param_name,'ss00')==1
        neuron(param);
    elseif strcmp(param_name,'ss01')==1
%		neuron_network(param);
		neuron_ss01(param);
    elseif strcmp(param_name,'pattern_positive')==1
		neuron_pattern_positive(param);
	elseif strcmp(param_name,'sparse_coding')==1
        if param.fast_sparse==1
            if param.inter_neuron==0
                neuron_sparse_coding_fast(param);
            else
                neuron_sparse_coding_inter(param);
            end
        else
            neuron_sparse_coding(param);
        end
    elseif strcmp(param_name,'ss00toy')==1
        neuron_toy(param);
    elseif strcmp(param_name,'pattern_negative')==1
		fprintf('I will finish it!');
		return;
    else
        fprintf('Wrong Argument\n');
        return;
    end
end
