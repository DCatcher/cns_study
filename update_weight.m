function param   = update_weight(param)

for i=1:param.neuron_n
    tmp_fire_list           = find(param.fire_list(i,:)==1);
    for j=1:param.sig_n
        tmp_sig_in_list     = find(param.sig_in_list(j,:)==1);
        
        for k=1:length(tmp_fire_list)
            tmp_fire_time                       = tmp_fire_list(k);
            tmp_sig_time_list                   = (tmp_sig_in_list >= tmp_fire_time) & (tmp_sig_in_list < (tmp_fire_time + length(param.stdp_sig_to_ex_neg)));
            param.delta_g_a_sig_to_ex(i,j)      = param.delta_g_a_sig_to_ex(i,j) + sum(param.stdp_sig_to_ex_neg(tmp_sig_in_list(tmp_sig_time_list) - tmp_fire_time + 1));
            tmp_sig_time_list                   = (tmp_sig_in_list < tmp_fire_time) & (tmp_sig_in_list > (tmp_fire_time - length(param.stdp_sig_to_ex_neg)));
            param.delta_g_a_sig_to_ex(i,j)      = param.delta_g_a_sig_to_ex(i,j) + sum(param.stdp_sig_to_ex_pos(tmp_fire_time - tmp_sig_in_list(tmp_sig_time_list)));
        end
    end
end

for i=1:param.neuron_n
    tmp_fire_list           = find(param.fire_list(i,:)==1);
    for j=1:param.neuron_n
        tmp_sig_in_list     = find(param.fire_list(j,:)==1);
        
        for k=1:length(tmp_fire_list)
            tmp_fire_time                       = tmp_fire_list(k);
            tmp_sig_time_list                   = (tmp_sig_in_list >= tmp_fire_time) & (tmp_sig_in_list < (tmp_fire_time + length(param.stdp_ex_to_ex_neg)));
            param.delta_g_a_ex_to_ex(i,j)       = param.delta_g_a_ex_to_ex(i,j) + sum(param.stdp_ex_to_ex_neg(tmp_sig_in_list(tmp_sig_time_list) - tmp_fire_time + 1));
            tmp_sig_time_list                   = (tmp_sig_in_list < tmp_fire_time) & (tmp_sig_in_list > (tmp_fire_time - length(param.stdp_ex_to_ex_neg)));
            param.delta_g_a_ex_to_ex(i,j)      = param.delta_g_a_ex_to_ex(i,j) + sum(param.stdp_ex_to_ex_pos(tmp_fire_time - tmp_sig_in_list(tmp_sig_time_list)));
        end
    end
end

if mod(param.time_now, param.ln_time)==0
    param.g_a_sig_to_ex         = param.g_a_sig_to_ex + param.delta_g_a_sig_to_ex;
    param.g_a_sig_to_ex         = max(0, param.g_a_sig_to_ex); 
    param.delta_g_a_sig_to_ex   = zeros(param.neuron_n, param.sig_n);
    
    param.g_a_ex_to_ex          = param.g_a_ex_to_ex + param.delta_g_a_ex_to_ex;
    param.g_a_ex_to_ex          = max(0, param.g_a_ex_to_ex); 
    param.delta_g_a_ex_to_ex    = zeros(param.neuron_n, param.neuron_n);
end