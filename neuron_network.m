function neuron_network(param)
	sig_n = param.sig_n;
	ex_n = param.ex_n;
	in_m = param.m;
	tmax = param.tmax;
	fs = param.fs;
	conn_rate = param.conn_rate;
    time_simu = param.time_simu;
    cont = param.cont;
    
    expand_size = param.expand_size;
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
    
    short_report_mode = param.short_report_mode;

    time_all = fs*tmax;
    vol_ex = zeros(ex_n,time_all);
	vol_ex(:,1) = v_reset;
    vol_in = zeros(in_m,time_all);
    vol_in(:,1) = v_reset;
	time_all = tmax*fs;
	for lamda_now = lamda
        tag_now = round(1/lamda_now);
		if (cont==0)
			g_a_sig_to_ex = g_max*rand(ex_n, sig_n);
			p_a_sig_to_ex = zeros(ex_n, sig_n);
			M_sig_to_ex = 0*ones(ex_n,1);
			g_sig_to_ex = 1*ones(ex_n,1);
			conn_sig_to_ex = (rand(ex_n,sig_n)<conn_rate);
			
			g_a_ex_to_in = g_max*rand(in_m, ex_n);
			p_a_ex_to_in = zeros(in_m, ex_n);
			M_ex_to_in = 0*ones(in_m,1);
			g_ex_to_in = 1*ones(in_m,1);
			conn_ex_to_in = (rand(in_m, ex_n)<conn_rate);
			
			g_a_in_to_other = g_max_in*rand(ex_n+in_m, in_m);
			p_a_in_to_other = zeros(ex_n+in_m, in_m);
			M_in_to_other = 0*ones(ex_n+in_m,1);
			g_in_to_other = 0*ones(ex_n+in_m,1);
			conn_in_to_other = (rand(ex_n+in_m, in_m)<conn_rate);
		elseif (cont==1)
			load save_data
		end
		for time_now=0:tmax:time_simu
            if short_report_mode==1
                fprintf('time_now:%i/%i, step: %i, Poi fre: %i\n',time_now,time_simu,tmax,tag_now);
            end
			[sig_ex_all,sig_in_gen] = gen_sig_ex_1(sig_n,in_m,tmax,fs,lamda_now);
			if (time_now>0)
				vol_ex(:,1) = vol_ex(:,time_all);
				vol_in(:,1) = vol_in(:,time_all);
			end
			for i=2:time_all
				vol_ex(:,i) = -60*(vol_ex(:,i-1)>v_th);
				vol_ex(:,i) = vol_ex(:,i) + (vol_ex(:,i-1)<v_th).*...
					((v_rest-vol_ex(:,i-1)+(g_sig_to_ex.*(e_ex-vol_ex(:,i-1))+g_in_to_other(1:ex_n).*(e_in-vol_ex(:,i-1))))/tao_m*1.0/fs+vol_ex(:,i-1));
				Po_ex = (vol_ex(:,i)>=v_th);
				vol_ex(:,i) = vol_ex(:,i).*(vol_ex(:,i)<v_th); 

%sig_to_ex begin
				sig_in = sig_ex_all(i,:);
				sig_in_all = repmat(sig_in, ex_n, 1);
				sig_in_all = conn_sig_to_ex.* sig_in_all;
				
				g_sig_to_ex = g_sig_to_ex+(-g_sig_to_ex)/tao_ex*1.0/fs;
				M_sig_to_ex = M_sig_to_ex + (-M_sig_to_ex)/tao_neg*1.0/fs;
				p_a_sig_to_ex = p_a_sig_to_ex+(-p_a_sig_to_ex)/tao_pos*1.0/fs;

				M_sig_to_ex = M_sig_to_ex - Po_ex*A_neg;
				g_sig_to_ex = g_sig_to_ex + (sum((g_a_sig_to_ex.*sig_in_all)'))';
				g_a_sig_to_ex = g_a_sig_to_ex + sig_in_all.*repmat(M_sig_to_ex, 1, sig_n)*g_max;
				p_a_sig_to_ex = p_a_sig_to_ex + A_pos*sig_in_all;

				g_a_sig_to_ex = g_a_sig_to_ex + p_a_sig_to_ex*g_max.*repmat(Po_ex, 1, sig_n);

				g_a_sig_to_ex = max(g_a_sig_to_ex,0);
				g_a_sig_to_ex = min(g_a_sig_to_ex,g_max);
%sig_to_ex finish

				vol_in(:,i) = -60*(vol_in(:,i-1)>v_th);
				tmp_vol_in = (vol_in(:,i-1)<v_th).*...
					((v_rest-vol_in(:,i-1)+(g_ex_to_in.*(e_ex-vol_in(:,i-1))+g_in_to_other((ex_n+1):(ex_n+in_m)).*(e_in-vol_in(:,i-1))))/tao_m*1.0/fs+vol_in(:,i-1));
				vol_in(:,i) = vol_in(:,i) + tmp_vol_in;
				Po_in = (vol_in(:,i)>=v_th);
				vol_in(:,i) = vol_in(:,i).*(vol_in(:,i)<v_th); 
%ex_to_in begin
				sig_in = vol_ex(:,i)'>v_th;
				sig_in_all = repmat(sig_in, in_m, 1);
				sig_in_all = conn_ex_to_in.* sig_in_all;
				
				g_ex_to_in = g_ex_to_in+(-g_ex_to_in)/tao_ex*1.0/fs;
				M_ex_to_in = M_ex_to_in + (-M_ex_to_in)/tao_neg*1.0/fs;
				p_a_ex_to_in = p_a_ex_to_in+(-p_a_ex_to_in)/tao_pos*1.0/fs;

				M_ex_to_in = M_ex_to_in - Po_in*A_neg;
				g_ex_to_in = g_ex_to_in + (sum((g_a_ex_to_in.*sig_in_all)'))';
				g_a_ex_to_in = g_a_ex_to_in + sig_in_all.*repmat(M_ex_to_in, 1, sig_n)*g_max;
				p_a_ex_to_in = p_a_ex_to_in + A_pos*sig_in_all;

				g_a_ex_to_in = g_a_ex_to_in + p_a_ex_to_in*g_max.*repmat(Po_in, 1, sig_n);

				g_a_ex_to_in = max(g_a_ex_to_in,0);
				g_a_ex_to_in = min(g_a_ex_to_in,g_max);                    
%ex_to_in finish                    

%in_to_other begin
				Po_all = [Po_ex;Po_in];

				sig_in = vol_in(:,i)'>v_th;
				sig_in_all = repmat(sig_in, ex_n+in_m, 1);
				sig_in_all = conn_in_to_other.* sig_in_all;
				
				g_in_to_other = g_in_to_other+(-g_in_to_other)/tao_ex_in*1.0/fs;
				M_in_to_other = M_in_to_other + (-M_in_to_other)/tao_neg_in*1.0/fs;
				p_a_in_to_other = p_a_in_to_other+(-p_a_in_to_other)/tao_pos_in*1.0/fs;

				M_in_to_other = M_in_to_other + Po_all*A_neg_in;
				g_in_to_other = g_in_to_other + (sum((g_a_in_to_other.*sig_in_all)'))';
				g_a_in_to_other = g_a_in_to_other + sig_in_all.*repmat(M_in_to_other-0.4*A_neg_in,1,in_m)*g_max_in;
				p_a_in_to_other = p_a_in_to_other + A_pos_in*sig_in_all;

				g_a_in_to_other = g_a_in_to_other + (p_a_in_to_other-0.4*A_pos_in)*g_max_in.*repmat(Po_all, 1, in_m);

				g_a_in_to_other = max(g_a_in_to_other,0);
				g_a_in_to_other = min(g_a_in_to_other,g_max_in);
%in_to_other finish

			end
			ans_my = mean(sum(([vol_ex;vol_in]>-10)));
            if short_report_mode==1
                fprintf('spiking_time:%f\n',ans_my);
            end
		end
		sig_ex_all = [];
		sig_in_gen = [];
        s = ['data_',int2str(tag_now)];
        save(s);
        fprintf('saved!%s\n',s);
	end
end
