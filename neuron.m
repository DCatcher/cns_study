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
            g_a = g_max*rand(1, n);
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
                fprintf('time_now:%i/%i, step: %i, Poi fre: %i\n',time_now,time_simu,tmax,tag_now);
            end
            [sig_ex_all,sig_in_all] = gen_sig_ex_1(n,m,tmax,fs,lamda_now);
            if (time_now>0)
                vol(1) = vol(time_all);
            end
            for i=2:time_all
                Po_ex = 0;
                if (vol(i-1)>v_th)
                    vol(i) = -60;
                else
                    vol(i) = (v_rest-vol(i-1)+(g_ex*(e_ex-vol(i-1))+g_in*(e_in-vol(i-1))))/tao_m*1.0/fs+vol(i-1);
                    if (vol(i)>=v_th)
                        Po_ex = 1;
                        vol(i) = 0;
                    end
                end

                sig_in = sig_ex_all(i,:);
                g_ex = g_ex+(-g_ex)/tao_ex*1.0/fs;
                M = M + (-M)/tao_neg*1.0/fs;
                p_a = p_a+(-p_a)/tao_pos*1.0/fs;

                M = M - Po_ex*A_neg;
                g_ex = g_ex + sum(g_a.*sig_in);
                g_a = g_a + sig_in*M*g_max;
                p_a = p_a + A_pos*sig_in;

                g_a = g_a + p_a*g_max*Po_ex;

                g_a = max(g_a,0);
                g_a = min(g_a,g_max);

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
            
            ans_my = sum(vol>-10);
            if short_report_mode==1
                fprintf('spiking_time:%i\n',ans_my);
            end

            if display_mode==1
                hist(g_a/g_max,M_s);
                title('g_a');
                pause(0.01);
            end
        end
        sig_ex_all = [];
        sig_in_all = [];
        s = ['data_',int2str(tag_now)];
        save(s);
        fprintf('saved!%s\n',s);
    end
end
