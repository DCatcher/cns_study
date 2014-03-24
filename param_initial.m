function param = param_initial(param_name)
    param.fs = 1000;
    param.display_mode = 1;
    param.short_report_mode = 1;
    
    if strcmp(param_name,'ss00')==1
        param.n = 1000;
        param.m = 200;
        param.tmax = 10;
        param.time_simu = 2000;
        
        param.tao_m = 20.0/1000;
        param.v_rest = -70;
        param.e_ex = 0;
        param.e_in = -70;
        param.v_th = -54;
        param.v_reset = -60;
        
        param.tao_ex = 5.0/1000;
        param.g_max = 0.015;
        param.tao_neg = 20.0/1000;
        param.tao_pos = 20.0/1000;
        param.A_pos = 0.005;
        param.A_neg = 1.05*param.A_pos;

        param.tao_ex_in = 5.0/1000;
        param.g_max_in = 0.5;
        param.tao_neg_in = 30.0/1000;
        param.tao_pos_in = 30.0/1000;
        param.A_pos_in = 0.015;
        param.A_neg_in = param.A_pos_in;        
        param.g_in_ba = 0.05;
        
        param.lamda = [1.0/10 1.0/15 1.0/20 1.0/25];
        
        param.cont = 0;
    elseif strcmp(param_name,'network')==1
		param.sig_n = 100;
		param.ex_n = 100;
		param.m = 20;
		param.tmax = 5;
		param.conn_rate = 0.8;
		param.time_simu = 2000;
		

		param.expand_size = 7;
		param.tao_m = 20.0/1000;
		param.v_rest = -70;
		param.e_ex = 0;
		param.e_in = -70;
		param.v_th = -54;
		param.v_reset = -60;
		
		param.tao_ex = 5.0/1000;
		param.g_max = 0.015*param.expand_size;
		param.tao_neg = 20.0/1000;
		param.tao_pos = 20.0/1000;
		param.A_pos = 0.005*param.expand_size;
		param.A_neg = 1.05*param.A_pos;
		
		param.tao_ex_in = 5.0/1000;
		param.g_max_in = 0.2*param.expand_size;
		param.tao_neg_in = 30.0/1000;
		param.tao_pos_in = 30.0/1000;
		param.A_pos_in = 0.015*param.expand_size;
		param.A_neg_in = param.A_pos_in;
        
		param.lamda = [1.0/10 1.0/15 1.0/20 1.0/25];
        param.cont = 0;
    elseif strcmp(param_name,'pattern_positive')==1
    elseif strcmp(param_name,'pattern_negative')==1
    else
    end
    
end
