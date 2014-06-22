function [vol, g_a] = neuron_toy(param)
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
	lamda_ex = param.lamda_ex;
    lamda_in = param.lamda_in;
    
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

	train_and_test = param.train_and_test;
    
    display_mode = param.display_mode;
    short_report_mode = param.short_report_mode;

	one_fire = param.one_fire;
	if one_fire==1
		fire_interval = param.fire_interval;
		fire_frame_len = param.fire_frame_len;
	end

	time_all = tmax*fs;
	vol = zeros(1,time_all);
	vol(1) = v_reset;
    M_s = 0:0.05:1;
                
    if display_mode==1
        figure_bar = figure;
        figure_plot = figure;
    end
    
	g_ave_all = [];
	rate_ave_all = [];

    for rate_A_np = 1.05:0.01:1.05
        A_neg = A_pos*rate_A_np;
        rate_ave = [];
		g_ave = [];
        
        for lamda_now = 20:3:20
		if train_and_test==1
			pt_fre = param.pattern_fre;
			len_pat = length(pt_fre);
			lamda_per_fre = floor(fs/len_pat);
			lamda_ex = zeros(1, n);
			for i=1:len_pat
				lamda_ex(((i-1)*lamda_per_fre+1):(i*lamda_per_fre)) = 1.0/pt_fre(i)*ones(1, lamda_per_fre);
			end

			lamda_in = zeros(1, m);
			for i=1:len_pat
				lamda_in(((i-1)*lamda_per_fre+1):(i*lamda_per_fre)) = 1.0/pt_fre(i)*ones(1, lamda_per_fre);
            end
            
            lamda_in = param.lamda_in*10/10;
		else
	        lamda_ex = param.lamda_ex*1.0/20;
			lamda_in = param.lamda_in*10/10;
        end
        
        if param.mode_pic==1
            sigma = 0.5;
            tao_c = 0.02;
            c_a = zeros(1, n);
            
			pt_fre = param.pattern_fre;
			len_pat = length(pt_fre);
			lamda_per_fre = floor(fs/len_pat);
            for i=1:len_pat
				c_a(((i-1)*lamda_per_fre+1):(i*lamda_per_fre)) = pt_fre(i)*ones(1, lamda_per_fre)/100*0.4;
            end 
            sigma_a = sqrt(sigma^2-c_a.^2);
        end
        spike_ave = [];
        if (cont==0)
            g_a = g_max*rand(1, n);
            p_a = zeros(1, n);
            M = 0;
            g_a_in = g_max_in*rand(1, m);
            p_a_in = zeros(1, m);
            M_in = 0;
            g_ex = 1;
            g_in = 0;
        elseif (cont==1)
            load save_data
        end
        tag_now = lamda_now;
        
        for time_now=0:tmax:time_simu
			g_ex = 0; 
			g_in = 0;

            if short_report_mode==1
                fprintf('time_now:%f/%f, step: %f, Poi fre: %f, Rate:%f\n',time_now,time_simu,tmax,tag_now, rate_A_np);
            end
            if param.mode_pic==1
                [sig_ex_all,sig_in_all] = gen_sig_ex_corr(n,m,tmax,fs, 1.0/20, c_a, tao_c, sigma_a);
            elseif one_fire==0
                [sig_ex_all,sig_in_all] = gen_sig_ex_toy(n,m,tmax,fs, lamda_ex, lamda_in);
			else
				[sig_ex_all,sig_in_all] = gen_sig_one_fire(n,m,tmax,fs,param.pattern_fre, fire_interval);
            end

            if (time_now>0)
				if one_fire==0
	                vol(1) = vol(time_all);
				else
					vol(1) = v_reset;
					g_ex = 0;
					M = 0;
            		p_a = zeros(1, n);
					g_in = 0;
            		p_a_in = zeros(1, m);
					M_in = 0;
				end
            end
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
% 	            g_a_in = g_a_in + sig_in*(M_in-param.A_neg_max*A_neg_in)*g_max_in;
%                 p_a_in = p_a_in + A_pos_in*sig_in;
% 
%                	g_a_in = g_a_in + (p_a_in-param.A_neg_max*A_pos_in)*g_max_in*Po_ex;
% 
%                 g_a_in = max(g_a_in,0);
%                 g_a_in = min(g_a_in,g_max_in);
            end
            
            ans_my = sum(vol>-10);
            spike_ave(end+1) = ans_my;
            if short_report_mode==1
                fprintf('spiking time:%i\n',ans_my);
            end
 
 			hist(g_a/g_max, M_s);
 			title('g_{ex}');
 			pause(0.01);

            if display_mode==1
				if param.mode_pic==0
                    figure(figure_bar);
                    g_left =g_a((length(param.tmp_lamda)+1):length(g_a)); 
 					hist(g_left/g_max,M_s);
                 	title('g_{left}');
% 					plot(vol);
                    figure(figure_plot);
                    g_part = g_a(1:length(param.tmp_lamda));
                    hist(g_part/g_max, M_s);
                    title('g_{part}');
                    fprintf('spiking weight rate:%f\n',mean(g_left)/mean(g_part));
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
%         sig_ex_all = [];
%         sig_in_all = [];
%         s = ['data_',int2str(tag_now)];
%         save(s);
%         fprintf('saved!%s\n',s);
        rate_ave(end+1) = mean(spike_ave(end-param.ave_range:end));
		g_ave(end+1) = mean(g_a/g_max);
        end
        %figure;
        %plot(rate_ave);
        %pause(0.01);

		g_ave_all = [g_ave_all; g_ave];
		rate_ave_all = [rate_ave_all; rate_ave];

		sta_pat = zeros(1, len_pat);
		sta_num = 0;

		mean_fre = mean(param.pattern_fre);
		sum_fre = sum(abs(param.pattern_fre-mean_fre));
		if train_and_test==1
			for test_num_now = 0:param.test_num
				pt_fre = 100*rand(1, len_pat);
				pt_fre = pt_fre-mean(pt_fre);
				pt_fre = pt_fre/sum(abs(pt_fre))*sum_fre;
				pt_fre = pt_fre+mean_fre;
				pt_fre = max(pt_fre, 1);
				pt_fre = min(pt_fre, 100);

				lamda_ex = zeros(1, n);
				for i=1:len_pat
					lamda_ex(((i-1)*lamda_per_fre+1):(i*lamda_per_fre)) = 1.0/pt_fre(i)*ones(1, lamda_per_fre);
                end

                if param.mode_pic==1
                    sigma = 0.5;
                    tao_c = 0.02;
                    c_a = zeros(1, n);

                    for i=1:len_pat
                        c_a(((i-1)*lamda_per_fre+1):(i*lamda_per_fre)) = pt_fre(i)*ones(1, lamda_per_fre)/100*0.4;
                    end 
                    sigma_a = sqrt(sigma^2-c_a.^2);
                end                

				lamda_in = zeros(1, m);
				for i=1:len_pat
					lamda_in(((i-1)*lamda_per_fre+1):(i*lamda_per_fre)) = 1.0/pt_fre(i)*ones(1, lamda_per_fre);
                end
                
                lamda_in = param.lamda_in*10/10;

				if short_report_mode==1
					fprintf('test num now:%i/%i, step: %i\n',test_num_now,param.test_num,1);
                end
                
                if param.mode_pic==1
                    [sig_ex_all,sig_in_all] = gen_sig_ex_corr(n,m,tmax,fs, 1.0/20, c_a, tao_c, sigma_a);
				elseif one_fire==0
                    [sig_ex_all,sig_in_all] = gen_sig_ex_toy(n,m,tmax,fs, lamda_ex, lamda_in);
				else
					[sig_ex_all,sig_in_all] = gen_sig_one_fire(n,m,tmax,fs,param.pattern_fre, fire_interval);
                end
                
				vol = zeros(1, time_all);
				vol(1) = v_reset;
				g_ex = 0;
				g_in = 0;
				M = 0;
				M_in = 0;
            	p_a = zeros(1, n);
            	p_a_in = zeros(1, m);

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
						end
					end

					sig_in = sig_ex_all(i,:);
					g_ex = g_ex+(-g_ex)/tao_ex*1.0/fs;
					g_ex = g_ex + sum(g_a.*sig_in);

					sig_in = sig_in_all(i,:);
                	g_in = g_in +(-g_in)/tao_ex_in*1.0/fs;
					g_in = g_in + sum(g_a_in.*sig_in);
				end
				spike_time = sum(vol>-10);
				sta_pat = sta_pat + spike_time*pt_fre;
				sta_num = sta_num + spike_time;
				plot(sta_pat/sta_num);
				hold on
				plot(param.pattern_fre,'r')
				hold off
				pause(0.01);
			end
		end
    end
    
    sig_ex_all = [];
    sig_in_all = [];
    s = 'data_toy_test';
    save(s);
    fprintf('saved!%s\n',s);
end
