function [vol, g_a] = neuron(param)
    n = param.n;
    m = param.m;
    tmax = param.tmax;
    fs = param.fs;
    time_simu = param.time_simu;
    cont = param.cont;

	tao_m = param.tao_m;
	v_rest = param.v_rest;
	e_ex = param.e_ex;
	e_in = param.e_in;
	v_th = param.v_th;
	v_reset = param.v_reset;
	lamda = param.lamda;
    
	tao_ex = param.tao_ex;
	g_max = param.g_max;
	tao_neg = param.tao_neg;
	tao_pos = param.tao_pos;
	A_pos = param.A_pos;
	A_neg = param.A_neg;
    
    tao_ex_in = param.tao_ex_in;
	g_max_in = param.g_max_in;
	tao_neg_in = param.tao_neg_in;
	tao_pos_in = param.tao_pos_in;
	A_pos_in = param.A_pos_in;
	A_neg_in = param.A_neg_in;
    g_in_ba = param.g_in_ba;
    
    display_mode = param.display_mode;
    short_report_mode = param.short_report_mode;

	time_all = tmax*fs;
	vol = zeros(1,time_all);
	vol(1) = v_reset;
                
    if display_mode==1
        M_s = 0:0.05:1;
        figure;
    end

    for lamda_now=lamda
        if (cont==0)
            if param.fix_g_a==0
                g_a = g_max*rand(1, n);
            else
                g_a = [g_max*param.init_g_a*ones(1, n/2), g_max*param.init_g_a_2*ones(1, n/2)];
            end
            p_a = zeros(1, n);
            M = 0;
            g_a_in = g_max*zeros(1, m);
            p_a_in = zeros(1, m);
            M_in = 0;
            g_ex = 1;
            g_in = 0;
        elseif (cont==1)
            load save_data
        end
        tag_now = round(1/lamda_now);
        
        for time_now=0:tmax:time_simu
            if short_report_mode==1
                fprintf('time_now:%i/%i, step: %i, Poi fre: %i, ',time_now,time_simu,tmax,tag_now);
            end
			if param.mode_pic==0
                if param.with_speed==0                    
                        [sig_ex_all,sig_in_all] = gen_sig_ex_1(n,m,tmax,fs,lamda_now);
                else
                    if param.sep_sig==0
                        [sig_ex_all,sig_in_all, speed_ex] = gen_sig_ex_with_speed(n,m,tmax,fs,lamda_now);
                    else
                        [sig_ex_all_1,sig_in_all, speed_ex_1] = gen_sig_ex_with_speed(n/2,m,tmax,fs,lamda_now);
                        [sig_ex_all_2,sig_in_all_tmp, speed_ex_2] = gen_sig_ex_with_speed(n/2,0,tmax,fs,param.an_lamda);
                        sig_ex_all  = [sig_ex_all_1(1:tmax*fs, :), sig_ex_all_2(1:tmax*fs, :)];
                        speed_ex    = [speed_ex_1(1:tmax*fs, :), speed_ex_2(1:tmax*fs, :)];
                    end
                end
			else
				[sig_ex_all,sig_in_all] = gen_sig_ex_corr(n,m,tmax,fs,lamda_now, param.c_a, param.tao_c, param.sigma_a);
            end
            if (time_now>0)
                vol(1) = vol(time_all);
            end
            
            if param.with_speed==1
                vol_sig     = v_reset*ones(1, n);
                speed_in    = zeros(1, n);
            end

            fire_list = [];
            sig_list = cell(n, 1);
            
            for i=2:time_all
                Po_ex = 0;
                if (vol(i-1)>v_th)
                    vol(i) = -60;
                else
                    vol(i) = (v_rest-vol(i-1)+(g_ex*(e_ex-vol(i-1))+g_in*(e_in-vol(i-1))))/tao_m*1.0/fs+vol(i-1);
%                    vol(i) = (v_rest-vol(i-1)+(g_ex*(e_ex-vol(i-1))+g_in*(e_in-vol(i-1))))/tao_m*1.0/fs+vol(i-1);
                    if (vol(i)>=v_th)
                        Po_ex = 1;
                        vol(i) = 0;
                        fire_list(end+1)    = i;      
                        if param.with_speed==1
                            speed_in    = speed_in + g_a/g_max*param.neg_speed;
                        end
                    end
                end

                if param.with_speed==1
                    vol_sig             = (speed_ex(i,:) - speed_in)/tao_m*1.0/fs+vol_sig;
                    sig_ex_all(i, :)    = (vol_sig>=v_th);
                    speed_in            = (-speed_in)/param.neg_speed_tao*1.0/fs+speed_in;
                    vol_sig(vol_sig>=v_th)  = v_reset;
                end
                sig_in = double(sig_ex_all(i,:));
                g_ex = g_ex+(-g_ex)/tao_ex*1.0/fs;
                M = M + (-M)/tao_neg*1.0/fs;
                p_a = p_a+(-p_a)/tao_pos*1.0/fs;

                M = M - Po_ex*A_neg;
                g_ex = g_ex + sum(g_a.*sig_in);
                p_a = p_a + A_pos*sig_in;
                if param.fix_g_a==0
                    g_a = g_a + sig_in*M*g_max;
                    g_a = g_a + p_a*g_max*Po_ex;
                    g_a = max(g_a,0);
                    if param.without_limit==0
                        g_a = min(g_a,g_max);
                    end
                end

                sig_in = sig_in_all(i,:);
                g_in = (-g_in)/tao_ex_in*1.0/fs+g_in;
                g_in = g_in + sum(g_in_ba.*sig_in);
%                 g_in = g_in+(-g_in)/tao_ex_in*1.0/fs;
%                 M_in = M_in + (-M_in)/tao_neg_in*1.0/fs;
%                 p_a_in = p_a_in+(-p_a_in)/tao_pos_in*1.0/fs;
% 
%                 M_in = M_in + Po_ex*A_neg_in;
%                 g_in = g_in + sum(g_a_in.*sig_in);
%                 g_a_in = g_a_in + sig_in*(M_in-0.4*A_neg_in)*g_max_in;
%                 p_a_in = p_a_in + A_pos_in*sig_in;
% 
%                 g_a_in = g_a_in + (p_a_in-0.4*A_pos_in)*g_max_in*Po_ex;
% 
%                 g_a_in = max(g_a_in,0);
%                 g_a_in = min(g_a_in,g_max_in);

            end
            
            for i=1:n
                sig_list{i}     = find(sig_ex_all(:,i));
            end            
            
            ans_my = sum(vol>-10);
            if short_report_mode==1
                fprintf('spiking_time:%i\n',ans_my);
            end

            if display_mode==1
				if param.mode_pic==0
                    if param.without_limit==0
                        hist(g_a/g_max,M_s);
                    else
                        hist(g_a/g_max);
                    end
                	title('g_a');
					%plot(vol);
				else
					bar_len = 25;
					a = 1:bar_len;
					b = zeros(1,bar_len);
					for i=1:bar_len
						b(i) = mean(g_a(((i-1)*n/bar_len+1):(i*n/bar_len)));
					end
					bar(a,b/g_max);
					title('Fig 5.a');
				end
                pause(0.01);
            end
        end
        sig_ex_all = [];
        sig_in_all = [];
        s = [param.save_prefix, 'data_',int2str(tag_now)];
        save(s);
        fprintf('saved!%s\n',s);
    end
end
