function neuron_network(sig_n, ex_n,in_m,tmax,conn_rate, fs,cont)
    expand_size = 7;
    
	tao_m = 20.0/1000;
	v_rest = -70;
	e_ex = 0;
	e_in = -70;
	v_th = -54;
	v_reset = -60;
	lamda = 1.0/10;
	time_all = tmax*fs;
    
	tao_ex = 5.0/1000;
	g_max = 0.015*expand_size;
	tao_neg = 20.0/1000;
	tao_pos = 20.0/1000;
	A_pos = 0.005*expand_size;
	A_neg = 1.05*A_pos;
    
    tao_ex_in = 5.0/1000;
	g_max_in = 0.2*expand_size;
	tao_neg_in = 30.0/1000;
	tao_pos_in = 30.0/1000;
	A_pos_in = 0.015*expand_size;
	A_neg_in = A_pos_in;

	vol_ex = zeros(ex_n,time_all);
	vol_ex(:,1) = v_reset;
    vol_in = zeros(in_m,time_all);
    vol_in(:,1) = v_reset;
    big_time = 8;
	my_time = 359;%一个my_time就是tmax秒
% 	g_a_sig_to_ex = g_max*ones(ex_n, sig_n);
	%{
	g_a(1:230) = zeros(1,230);
	g_a(231:311) = 0.1*0.015*ones(1,81);
	g_a(312:411) = 0.9*0.015*ones(1,100);
	g_a(412:461) = 0.8*0.015*ones(1,50);
	%}
% 	if (cont==0)
% 		g_a = g_max*ones(1, n);
% 		p_a = zeros(1, n);
% 		M = 0;
% 		g_ex = 1;
% 		g_in = 0;
% 	elseif (cont==1)
% 		load save_data
% 	end

% 	g_ex_all=zeros(1,time_all);
% 	g_in_all=zeros(1,time_all);
% 	m_all=zeros(1,time_all);
% 	p_a_mean_all = zeros(1,time_all);

	%{
	if (cont==0)
		[sig_ex_all sig_in_all] = gen_sig_ex_1(n,m,tmax,fs,lamda);
		fprintf('Generate signal finished\n');
		save sig_1 sig_ex_all sig_in_all
	elseif (cont==1)
		load sig_1
	end
	%}
    for big_rand = 0:2
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
%         g_ex_all=zeros(1,time_all);
%         g_in_all=zeros(1,time_all);
%         m_all=zeros(1,time_all);
%         p_a_mean_all = zeros(1,time_all);
        
        if big_rand==0
            lamda = 1.0/15;
        elseif big_rand==1
            lamda = 1.0/20;
        elseif big_rand==2
            lamda = 1.0/25;
        end
        for big_time_now =0:big_time
            for time_now=0:my_time
				fprintf('time_now:%i/%i,round:%i/%i\n',time_now,my_time,big_time_now,big_time);
                [sig_ex_all,sig_in_gen] = gen_sig_ex_1(sig_n,in_m,tmax,fs,lamda);
                if (time_now>0)
                    vol_ex(:,1) = vol_ex(:,time_all);
                    vol_in(:,1) = vol_in(:,time_all);
                end
                for i=2:time_all
%                     disp(i);
%                     g_ex_all(i) = g_sig_to_ex;
%                     g_in_all(i) = g_in_to_other;
%                     m_all(i) = M_sig_to_ex;
%                     p_a_mean_all(i) = max(p_a);

%                     if (vol_ex(i-1)>v_th)
%                         vol_ex(i) = -60;
%                     else
%                         vol_ex(i) = (v_rest-vol_ex(i-1)+(g_sig_to_ex*(e_ex-vol_ex(i-1))+g_in_to_other*(e_in-vol_ex(i-1))))/tao_m*1.0/fs+vol_ex(i-1);
%                         if (vol_ex(i)>=v_th)
%                             Po_ex = 1;
%                             vol_ex(i) = 0;
%                         end
%                     end
                    
                    vol_ex(:,i) = -60*(vol_ex(:,i-1)>v_th);
                    vol_ex(:,i) = vol_ex(:,i) + (vol_ex(:,i-1)<v_th).*...
                        ((v_rest-vol_ex(:,i-1)+(g_sig_to_ex.*(e_ex-vol_ex(:,i-1))+g_in_to_other(1:ex_n).*(e_in-vol_ex(:,i-1))))/tao_m*1.0/fs+vol_ex(:,i-1));
                    Po_ex = (vol_ex(:,i)>=v_th);
                    vol_ex(:,i) = vol_ex(:,i).*(vol_ex(:,i)<v_th); 

%                     g_in = in_synp(g_in, sig_in_all(i,:), m, fs);
%                     [g_ex, M, p_a, g_a] = ex_synp(g_ex, M, p_a, g_a, Po_ex, sig_ex_all(i,:), n, fs);

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
%                     [g_in, M_in, p_a_in, g_a_in] = in_synp(g_in, M_in, p_a_in, g_a_in, Po_ex, sig_in_all(i,:), m, fs);

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
%                     sig_in = sig_in_gen(i,:);
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
    
                end
                ans_my = sum((vol_ex>-10));
 				fprintf('spiking_time:%f\n',mean(ans_my));
%                 plot(ans_my);
            end
            sig_ex_all = [];
            sig_in_gen = [];
			if big_rand==0
			 s = ['save_data_15_',int2str(big_time_now)];
			elseif big_rand==1
			 s = ['save_data_20_',int2str(big_time_now)];
			elseif big_rand==2
			 s = ['save_data_25_',int2str(big_time_now)];
			end
			save(s);
			fprintf('saved!%s\n',s);
    %         save [save_data] g_ex g_in M p_a g_a
        end
    end
	%{
	plot(g_ex_all);
	title('g_{ex}');
	figure;
	plot(p_a);
	title('p_{a}');
	figure;
	%}
end
