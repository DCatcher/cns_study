function cns(param_name, varargin)
    close all
    %rand('state',sum(100*clock));
    param   = param_initial(param_name);
    param   = parseArgs(varargin, param);
    param   = param_again(param);
    param   = parseArgs(varargin, param);
    
    rand('state',param.rand_seed);
	randn('state',param.randn_seed);
%     load IMAGES
%     param.IMAGES    = IMAGES;    
    
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
    elseif strcmp(param_name,'new')==1
		neuron_new(param);
    elseif strcmp(param_name,'pattern_negative')==1
		fprintf('I will finish it!');
		return;
    else
        fprintf('Wrong Argument\n');
        return;
    end
end
