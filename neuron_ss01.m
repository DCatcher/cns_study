function neuron_ss01(param)
	sig_n = param.sig_n;
	ex_n = param.ex_n;
	tmax = param.tmax;
	fs = param.fs;
    time_simu = param.time_simu;
    cont = param.cont;
	display_mode = param.display_mode;
    
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
	A_neg_recur = param.A_neg_recur;
    
    short_report_mode = param.short_report_mode;
	se_conn_rate = param.sig_to_ex_conn_rate;
	net_conn_rate = param.net_conn_rate;

    time_all = fs*tmax;
    vol_ex = zeros(ex_n,time_all);
	vol_ex(:,1) = v_reset;
	time_all = tmax*fs;

	for lamda_now = lamda
        tag_now = round(1/lamda_now);
		if (cont==0)
			for i=1:ex_n
				for j=1:sig_n
					d = j/(sig_n/ex_n)-i;
					if (d>100)
						d = 200-d;
					end
					if (d<-100)
						d = 200+d;
					end
					g_a_sig_to_ex(i,j) = 0.5*g_max*exp(-0.5*(d/100)^2);
				end
			end
			%g_a_sig_to_ex = rand(ex_n,sig_n)*g_max;
			%g_a_sig_to_ex(80:120,400:600) = rand(41,201)*g_max;
			conn_sig_to_ex = (rand(ex_n,sig_n)<se_conn_rate);
			p_a_sig_to_ex = zeros(ex_n, sig_n);
			M_sig_to_ex = 0*ones(ex_n,1);
			g_sig_to_ex = 1*ones(ex_n,1);
			%to be finished
			conn_net = max((rand(ex_n,ex_n)<net_conn_rate)-diag(ones(1,ex_n)),0);
			p_a_ex_to_ex = zeros(ex_n, ex_n);
			M_ex_to_ex = 0*ones(ex_n,1);
			g_ex_to_ex = 1*ones(ex_n,1);
			g_a_ex_to_ex = rand(ex_n,ex_n)*0.2*g_max;
		elseif (cont==1)
			load save_data
		end
		sig_back = repmat([0,1],1,5000);
		for time_now=0:tmax:time_simu
            if display_mode==1
				M_s = 30;
				%hist(g_a_sig_to_ex(1,:)/g_max,M_s);
				%{
				g_a = g_a_ex_to_ex(50,:);
					bar_len = 10;
					a = 1:bar_len;
					b = zeros(1,bar_len);
					for i=1:bar_len
						b(i) = mean(g_a(((i-1)*ex_n/bar_len+1):(i*ex_n/bar_len)));
					end
					bar(a,b/g_max);
				%}
%				subplot(2,1,1);
				g_a = 256-round(g_a_sig_to_ex/g_max*256);
				imshow(g_a,[1 256]);
				title('g_{se}');
				pause(0.02);
				%{
				subplot(2,1,2);
				g_a = 256-round(g_a_ex_to_ex/g_max*256);
				imshow(g_a,[1 256]);
				title('g_{ee}');
				pause(0.02);
				%}
			end
			
            if short_report_mode==1
				a_clk = clock;
                fprintf('%i:%i,%i/%i, step: %i, Poi fre: %i\n',a_clk(4),a_clk(5),time_now,time_simu,tmax,tag_now);
            end
			[sig_ex_all,sig_in_gen] = gen_sig_touch(sig_n,0,tmax,fs,lamda_now,param.r_1,param.tao_c,param.sigma);

			if (time_now>0)
				vol_ex(:,1) = vol_ex(:,time_all);
			end
			for i=2:time_all
				if (short_report_mode==1)
					if mod(i,round(time_all/20))==0
						fprintf('now %i/20\n',i/(round(time_all/20)));
					end
				end

				Po_ex = (vol_ex(:,i-1)>=v_th);
				vol_ex(:,i) = -60*(vol_ex(:,i-1)>v_th);
				vol_ex(:,i) = vol_ex(:,i) + (vol_ex(:,i-1)<v_th).*...
					((v_rest-vol_ex(:,i-1)+((g_sig_to_ex+g_ex_to_ex*(1-param.inhibi)+sig_back(i)*param.g_back).*(e_ex-vol_ex(:,i-1))))/tao_m*1.0/fs+param.g_in_stand*sum(Po_ex)*(param.inhibi)*(param.e_in-vol_ex(:,i-1))/tao_m*1.0/fs+vol_ex(:,i-1));
				Po_ex = (vol_ex(:,i)>=v_th);
				Po_ex_row = find(Po_ex==1);
				vol_ex(:,i) = vol_ex(:,i).*(vol_ex(:,i)<v_th); 

%sig_to_ex begin
				sig_in = sig_ex_all(i,:);
				sig_in_row = find(sig_in==1);
				
				g_sig_to_ex = g_sig_to_ex*(1-1/tao_ex*1.0/fs);
				M_sig_to_ex = M_sig_to_ex*(1-1/tao_neg*1.0/fs);
				p_a_sig_to_ex = p_a_sig_to_ex*(1-1/tao_pos*1.0/fs);

				M_sig_to_ex = M_sig_to_ex - Po_ex*A_neg;
				g_sig_to_ex = g_sig_to_ex + (sum((conn_sig_to_ex(:,sig_in_row).*g_a_sig_to_ex(:,sig_in_row))'))';
				g_a_sig_to_ex(:,sig_in_row) = g_a_sig_to_ex(:,sig_in_row) + repmat(M_sig_to_ex, 1, length(sig_in_row))*g_max;
				g_a_sig_to_ex(:,sig_in_row) = max(g_a_sig_to_ex(:,sig_in_row),0);
				p_a_sig_to_ex(:,sig_in_row) = p_a_sig_to_ex(:,sig_in_row) + A_pos;

				if length(Po_ex_row)<0.5*ex_n
					g_a_sig_to_ex(Po_ex_row,:) = g_a_sig_to_ex(Po_ex_row,:) + p_a_sig_to_ex(Po_ex_row,:)*g_max;
					g_a_sig_to_ex(Po_ex_row,:) = min(g_a_sig_to_ex(Po_ex_row,:),g_max);
				else
					g_a_sig_to_ex = g_a_sig_to_ex+p_a_sig_to_ex*g_max;
					Po_ex_row = find(Po_ex==0);
					g_a_sig_to_ex(Po_ex_row,:) = g_a_sig_to_ex(Po_ex_row,:) - p_a_sig_to_ex(Po_ex_row,:)*g_max;
					g_a_sig_to_ex = min(g_a_sig_to_ex,g_max);
					Po_ex_row = find(Po_ex==1);
				end

%sig_to_ex finish

%ex_to_ex begins
				g_ex_to_ex = g_ex_to_ex*(1-1/tao_ex*1.0/fs);
				M_ex_to_ex = M_ex_to_ex*(1-1/tao_neg*1.0/fs);
				p_a_ex_to_ex = p_a_ex_to_ex*(1-1/tao_pos*1.0/fs);

				M_ex_to_ex = M_ex_to_ex - Po_ex*A_neg_recur;
				g_ex_to_ex = g_ex_to_ex + (sum((conn_net(:,Po_ex_row).*g_a_ex_to_ex(:,Po_ex_row))'))';
				g_a_ex_to_ex(:,Po_ex_row) = g_a_ex_to_ex(:,Po_ex_row) + repmat(M_ex_to_ex, 1, length(Po_ex_row))*g_max;
				g_a_ex_to_ex(:,Po_ex_row) = max(g_a_ex_to_ex(:,Po_ex_row),0);
				p_a_ex_to_ex(:,Po_ex_row) = p_a_ex_to_ex(:,Po_ex_row) + A_pos;

				if length(Po_ex_row)<0.5*ex_n
					g_a_ex_to_ex(Po_ex_row,:) = g_a_ex_to_ex(Po_ex_row,:) + p_a_ex_to_ex(Po_ex_row,:)*g_max;
					g_a_ex_to_ex(Po_ex_row,:) = min(g_a_ex_to_ex(Po_ex_row,:),g_max);
				else
					g_a_ex_to_ex = g_a_ex_to_ex+p_a_ex_to_ex*g_max;
					Po_ex_row = find(Po_ex==0);
					g_a_ex_to_ex(Po_ex_row,:) = g_a_ex_to_ex(Po_ex_row,:) - p_a_ex_to_ex(Po_ex_row,:)*g_max;
					g_a_ex_to_ex = min(g_a_ex_to_ex,g_max);
					Po_ex_row = find(Po_ex==1);
				end
%ex_to_ex finishs
			end

            if short_report_mode==1
				a_clk = clock;
				fre_test = sum(sum((vol_ex>-10)));
                fprintf('%i:%i, Fre test: %i\n',a_clk(4),a_clk(5),fre_test);
            end
		end
		sig_ex_all = [];
		sig_in_gen = [];
        s = ['data_',int2str(tag_now)];
        save(s);
        fprintf('saved!%s\n',s);
	end
end
