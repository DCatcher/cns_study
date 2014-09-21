function param = param_again(param)

if strcmp(param.param_name, 'ss00')==1
end

if strcmp(param.param_name,'ss01')==1
    if param.input_2d==1
        param.sig_n = 1600;
        param.ex_n = 1;
        param.sigma = 5*sqrt(2);
    end
    
    if param.input_bar_2d==1
        param.gamma     = 3;
        param.sigma     = 8;
        param.sig_n     = 1600;
        param.ex_n      = 1;
    end    
    
    if param.input_bar_pic==1
        load IMAGES
        param.IMAGES    = IMAGES;
    end
end

if strcmp(param.param_name,'sparse_coding')==1

    if param.one_fire==1
        param.time_per = param.fire_interval*(param.sig_dim)^2+300;
        param.g_max = 0.45;
    end
    
        
    if param.one_fire_exp==1
        param.time_per = (max(param.fire_base_time_ex, param.fire_base_time_in) + param.fire_max_time) + 40;
        param.g_max = 0.1085;
%             param.g_max = 0.0685;
        param.offset = 2;
    end    
    
        
    if param.inter_neuron==1
        param.in_n = 10;
        param.g_max = 0.285;
        param.g_max_in = 100/(param.ex_n)*20*0.2*param.expand_size;
    end    
end

if strcmp(param.param_name,'ss00toy')==1

    if param.one_fire==1
        param.n = 100;
        param.fire_frame_len = param.fire_interval*param.n+200;
    else
        param.n = 1000;
    end    
    
    if param.one_fire==1
        param.tmax = param.fire_frame_len/param.fs;
    else     
        param.tmax = 10;
    end
    if param.one_fire==1
        param.time_simu = 200;
    else
        param.time_simu = 2000;
    end    
    
    if param.one_fire==1
        param.g_max = 0.35;
        param.A_pos = 0.05;
    else
        param.g_max = 0.015;
        param.A_pos = 0.005;
    end
    
    if param.one_fire==1
        param.pattern_fre = 100*rand(1, param.n);
    else
        param.pattern_fre = 100*rand(1, 20);
    end
    
end