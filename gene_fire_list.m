function param   = gene_fire_list(param)

param.sig_in_list   = zeros(param.sig_n, param.time_all);
param.fire_list     = zeros(param.neuron_n, param.time_all);

param.vol_ex        = param.v_reset * ones(param.neuron_n, param.time_all);

if param.with_speed==1
    param.vol_sig       = param.v_reset*ones(1, param.sig_n);
    param.speed_in      = zeros(1, param.sig_n);
end

for i=2:param.time_all
    if param.with_speed==1
        param.vol_sig           = (param.speed_ex(i,:) - param.speed_in)/param.tao_m*1.0/param.fs + param.vol_sig;
        sig_ex_all              = double(param.vol_sig>=param.v_th);
        param.speed_in          = (-param.speed_in)/param.neg_speed_tao*1.0/param.fs + param.speed_in;
        param.vol_sig( param.vol_sig >= param.v_th )    = param.v_reset;
    end
    
    sig_in_row          = sig_ex_all==1;
    param.sig_in_list(:,i)      = sig_in_row';
    
    param.vol_ex(:,i)   = param.v_reset * (param.vol_ex(:,i-1) > param.v_th);
    param.vol_ex(:,i)   = param.vol_ex(:,i) + (param.vol_ex(:,i-1) < param.v_th).* ...
        ((param.v_rest-param.vol_ex(:,i-1) + ...
        (param.g_sig_to_ex.*(param.e_ex - param.vol_ex(:,i-1))) + ...
        (param.g_ex_to_ex.*(param.e_in - param.vol_ex(:,i-1))))/param.tao_m*1.0/param.fs + ...
        param.vol_ex(:,i-1));
    Po_ex                   = (param.vol_ex(:,i) >= param.v_th);
    param.fire_list(:,i)    = Po_ex;
    
    param.vol_ex(:,i)   = param.vol_ex(:,i) .* (param.vol_ex(:,i) < param.v_th);
    param.vol_ex(:,i)   = max(param.e_in, param.vol_ex(:,i));
    
    param.g_sig_to_ex   = param.g_sig_to_ex * (1-1/param.tao_ex*1.0/param.fs);
    param.g_sig_to_ex   = param.g_sig_to_ex + (sum((param.g_a_sig_to_ex(:,sig_in_row))',1))';
    
    param.g_ex_to_ex    = param.g_ex_to_ex * (1-1/param.tao_ex*1.0/param.fs);
    param.g_ex_to_ex    = param.g_ex_to_ex + (sum((param.g_a_ex_to_ex(:,Po_ex))',1))';
    
    if param.with_speed==1
        for j=1:param.neuron_n
            if param.vol_ex(j, i) > param.v_th
                param.speed_in      = param.speed_in + param.g_a_sig_to_ex(j, :)/param.g_max * param.neg_speed;
            end
        end
    end
end