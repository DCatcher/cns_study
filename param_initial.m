function param = param_initial(param_name)
    param.fs = 1000;
    param.display_mode = 0;
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
        param.tao_neg = 32.7/1000;
        param.tao_pos = 16.8/1000;
        param.A_pos = 0.005;
        param.A_neg = 0.516*param.A_pos;

        param.tao_ex_in = 5.0/1000;
        param.g_max_in = 0.5;
        param.tao_neg_in = 30.0/1000;
        param.tao_pos_in = 30.0/1000;
        param.A_pos_in = 0.015;
        param.A_neg_in = param.A_pos_in;        
        param.g_in_ba = 0.05;
        
%        param.lamda = [1.0/10 1.0/15 1.0/20 1.0/25];
		param.lamda = [1.0/64];
        
        param.cont = 0;

		param.mode_pic = 0;

		param.sigma = 0.5;
		param.tao_c = 0.2;
		param.c_a = 0:0.0002:0.1998;
		param.sigma_a = sqrt(param.sigma^2-param.c_a.^2);
    elseif strcmp(param_name,'ss01')==1
		param.sig_n = 1000;
		param.ex_n = 200;
		param.m = 0;
		param.tmax = 10;
		param.sig_to_ex_conn_rate = 0.3;
		param.net_conn_rate = 0.5;
		param.time_simu = 2000;

		param.expand_size = 1;
		param.tao_m = 20.0/1000;
		param.v_rest = -74;
		param.e_ex = 0;
		param.e_in = -70;
		param.v_th = -54;
		param.v_reset = -60;
		
		param.tao_ex = 5.0/1000;
		param.g_max = 0.02;
		param.tao_neg = 20.0/1000;
		param.tao_pos = 20.0/1000;
		param.A_pos = 0.001*param.expand_size;
		param.A_neg = 1.06*param.A_pos;
		param.A_neg_recur = 1.04*param.A_pos;

		param.g_in_stand = 0.3*param.g_max;

		param.r_1 = 80;
		param.sigma = 100;
		param.tao_c = 20.0/1000;
		param.fre_back = 500;
		param.g_back = 0.096;
		
		%{
		param.tao_ex_in = 5.0/1000;
		param.g_max_in = 0.2*param.expand_size;
		param.tao_neg_in = 30.0/1000;
		param.tao_pos_in = 30.0/1000;
		param.A_pos_in = 0.015*param.expand_size;
		param.A_neg_in = param.A_pos_in;
		%}
        
%		param.lamda = [1.0/10 1.0/15 1.0/20 1.0/25];
		param.lamda = [1.0/10];
        param.cont = 0;
		param.inhibi = 0;
		param.ex_dist = 40;
    elseif strcmp(param_name,'pattern_positive')==1
        param.n = 1000;
        param.m = 200;
        param.tmax = 10;
        param.time_simu = 20000;
        
        param.tao_m = 20.0/1000;
        param.v_rest = -70;
        param.e_ex = 0;
        param.e_in = -70;
        param.v_th = -54;
        param.v_reset = -60;
        
        param.tao_ex = 5.0/1000;
        param.g_max = 0.015;
        param.tao_neg = 32.0/1000;
        param.tao_pos = 16.0/1000;
        param.A_pos = 0.005;
        param.A_neg = 0.525*param.A_pos;

        param.tao_ex_in = 5.0/1000;
        param.g_max_in = 0.5;
        param.tao_neg_in = 30.0/1000;
        param.tao_pos_in = 30.0/1000;
        param.A_pos_in = 0.015;
        param.A_neg_in = param.A_pos_in;        
        param.g_in_ba = 0.05;
        
        param.lamda = [1.0/35];

        param.cont = 0;

		param.time_pa = 0.05;
		param.wid_pa = 500;
		param.lamda_pa = 0.14;
    elseif strcmp(param_name,'sparse_coding')==1
		%note: threshold change?
		param.sig_dim = 8;
		param.ex_n = 64;
		param.sig_to_ex_conn_rate = 1;
		param.net_conn_rate = 1;
		param.train_num = 200000;
		param.tr_mat = 'IMAGES';
		param.lamda_max = 100;
        param.lamda_min = 10;
        param.cont = 0;
		param.time_per = 500;
		param.batch_size = 200;
        param.batch_learn = 1;
        param.easy_inhib = 0;
		param.time_pattern = 0;
        param.sta_big_num = 400;

		param.expand_size = 1;
		param.tao_m = 20.0/1000;
		param.v_rest = -74;
		param.e_ex = 0;
		param.e_in = -70;
		param.v_th = -54;
		param.v_reset = -60;
		
		param.tao_ex = 5.0/1000;
		param.g_max = 0.03;

		param.tao_neg = 20.0/1000;
		param.tao_pos = 20.0/1000;
		param.A_pos = 0.001*param.expand_size;
		%param.A_neg = 1.027*param.A_pos; %easy inhibi
		param.A_neg = 1.04*param.A_pos;
        
%         param.tao_neg = 32.0/1000;
%         param.tao_pos = 16.0/1000;
%         param.A_pos = 0.005;
%         param.A_neg = 0.51*param.A_pos;

		%param.tao_ex_in = 5.0/1000;
		param.tao_ex_in = 10.0/1000;
		param.g_max_in = 0.2*param.expand_size;
		param.tao_neg_in = 20.0/1000;
		param.tao_pos_in = 20.0/1000;
		param.A_pos_in = 0.005*param.expand_size;
		param.A_neg_in = param.A_pos_in;
		param.neg_amp = 0.3;
        
        param.inhibi_strength = 1*param.g_max_in;

		param.fre_min = 100;
		param.fre_max = 200;
		param.divide_len = 50;
		param.jitter_stride = 6;
        
%		param.lamda = [1.0/10 1.0/15 1.0/20 1.0/25];
    elseif strcmp(param_name,'ss00toy')==1
        param.n = 1000;
        param.m = 200;
        param.tmax = 10;
        param.time_simu = 2000;
        param.ave_range = 100;
        
        param.tao_m = 20.0/1000;
        param.v_rest = -70;
        param.e_ex = 0;
        param.e_in = -70;
        param.v_th = -54;
        param.v_reset = -60;
        
        param.tao_ex = 5.0/1000;
        param.g_max = 0.015;
        param.tao_neg = 20/1000;
        param.tao_pos = 20/1000;
        param.A_pos = 0.005;
        param.A_neg = 1.01*param.A_pos;

        param.tao_ex_in = 5.0/1000;
        param.g_max_in = 0.5;
        param.tao_neg_in = 30.0/1000;
        param.tao_pos_in = 30.0/1000;
        param.A_pos_in = 0.015;
        param.A_neg_in = param.A_pos_in;        
        param.g_in_ba = 0.05;
        
%        param.lamda = [1.0/10 1.0/15 1.0/20 1.0/25];
        param.tmp_lamda = 1*ones(1,200);
        base_fre = 1;
		param.lamda_ex = [1.0./param.tmp_lamda, 1.0/base_fre*ones(1, param.n-length(param.tmp_lamda))];
        param.lamda_in = 1.0/10*ones(1,param.m);
        
        param.cont = 0;

		param.mode_pic = 0;

		param.sigma = 0.5;
		param.tao_c = 0.2;
		param.c_a = 0:0.0002:0.1998;
		param.sigma_a = sqrt(param.sigma^2-param.c_a.^2);        
    elseif strcmp(param_name,'pattern_negative')==1
    else
    end
    
end
