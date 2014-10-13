function param   = short_report(param)

if param.short_report_mode==0
    return;
end

a_clk       = clock;
tmp_fire    = sum(sum(param.fire_list));
tmp_sig_in  = sum(sum(param.sig_in_list));
fprintf('%i:%i, now time: %.1f, fire fre: %.4f, sig fre: %.4f\n',...
    a_clk(4),a_clk(5),...
    param.time_now, tmp_fire/(param.neuron_n * param.tmax), tmp_sig_in/(param.sig_n * param.tmax));