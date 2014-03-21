function [vol, g_a] = neuron(n,m,tmax,fs,cont)
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
	A_pos = 0.005;
	A_neg = 1.05*A_pos;
    
    tao_ex_in = 5.0/1000;
	g_max_in = 0.2*expand_size;
	tao_neg_in = 30.0/1000;
	tao_pos_in = 30.0/1000;
	A_pos_in = 0.015;
	A_neg_in = A_pos_in;

	vol = zeros(1,time_all);
	vol(1) = v_reset;
    big_time = 8;
	my_time = 359;%一个my_time就是tmax秒
	g_a = g_max*ones(1, n);
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
            g_a = g_max*rand(1, n);
            p_a = zeros(1, n);
            M = 0;
            g_a_in = g_max*rand(1, m);
            p_a_in = zeros(1, m);
            M_in = 0;
            g_ex = 1;
            g_in = 0;
        elseif (cont==1)
            load save_data
        end
        g_ex_all=zeros(1,time_all);
        g_in_all=zeros(1,time_all);
        m_all=zeros(1,time_all);
        p_a_mean_all = zeros(1,time_all);
        
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
                [sig_ex_all,sig_in_all] = gen_sig_ex_1(n,m,tmax,fs,lamda);
                if (time_now>0)
                    vol(1) = vol(time_all);
                end
                for i=2:time_all
                    Po_ex = 0;
                    g_ex_all(i) = g_ex;
                    g_in_all(i) = g_in;
                    m_all(i) = M;
                    p_a_mean_all(i) = max(p_a);
                    if (vol(i-1)>v_th)
                        vol(i) = -60;
                    else
                        vol(i) = (v_rest-vol(i-1)+(g_ex*(e_ex-vol(i-1))+g_in*(e_in-vol(i-1))))/tao_m*1.0/fs+vol(i-1);
                        if (vol(i)>=v_th)
                            Po_ex = 1;
                            vol(i) = 0;
                        end
                    end
%                     g_in = in_synp(g_in, sig_in_all(i,:), m, fs);
%                     [g_ex, M, p_a, g_a] = ex_synp(g_ex, M, p_a, g_a, Po_ex, sig_ex_all(i,:), n, fs);

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
    
%                     [g_in, M_in, p_a_in, g_a_in] = in_synp(g_in, M_in, p_a_in, g_a_in, Po_ex, sig_in_all(i,:), m, fs);

                    sig_in = sig_in_all(i,:);
                    g_in = g_in+(-g_in)/tao_ex_in*1.0/fs;
                    M_in = M_in + (-M_in)/tao_neg_in*1.0/fs;
                    p_a_in = p_a_in+(-p_a_in)/tao_pos_in*1.0/fs;

                    M_in = M_in + Po_ex*A_neg_in;
                    g_in = g_in + sum(g_a_in.*sig_in);
                    g_a_in = g_a_in + sig_in*(M_in-0.4*A_neg_in)*g_max_in;
                    p_a_in = p_a_in + A_pos_in*sig_in;

                    g_a_in = g_a_in + (p_a_in-0.4*A_pos_in)*g_max_in*Po_ex;

                    g_a_in = max(g_a_in,0);
                    g_a_in = min(g_a_in,g_max_in);
    
                end
                ans_my = sum(vol>-10);
				fprintf('spiking_time:%i\n',ans_my);
            end
            sig_ex_all = [];
            sig_in_all = [];
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
