function neuron_new(param)

param.time_now      = 0;

disp(param);

while param.time_now < param.time_simu*param.fs
    
    param.time_now  = param.tmax_m + param.time_now;
    
    param   = gene_sig(param);
    param   = gene_fire_list(param);
    param   = update_weight(param);
    param   = display_new(param);
    param   = write_file(param);
    param   = short_report(param);
end