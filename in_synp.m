% function [g_in] = in_synp(g_in_o, sig_in, m, fs);
% %g_in_o是上一个时间点的值，sig_in是所有的in-synp的输入信号，m是in-synp的个数
% %fs是抽样频率（也就是模拟的最小时间的倒数）
% 
% 	tao_in = 5.0/1000;
% 	g_in_ba = 0.05;
% 
% 	g_in = (-g_in_o)/tao_in*1.0/fs+g_in_o;
% 
% 	for i=1:m
% 		if (sig_in(i)==1)
% 			g_in = g_in + g_in_ba;
% 		end
% 	end
% end

function [g_ex, M, p_a, g_a] = in_synp(g_ex_o, M_o, p_a_o, g_a_o, Po_ex, sig_in, n, fs);

	tao_ex = 5.0/1000;
	g_max = 0.015;
	tao_neg = 20.0/1000;
	tao_pos = 20.0/1000;
	A_pos = 0.005;
	A_neg = A_pos;

	g_ex = g_ex_o+(-g_ex_o)/tao_ex*1.0/fs;
	M = M_o + (-M_o)/tao_neg*1.0/fs;
	p_a = p_a_o+(-p_a_o)/tao_pos*1.0/fs;
	g_a = g_a_o;
	
    M = M + Po_ex*A_neg;
    g_ex = g_ex + sum(g_a.*sig_in);
    g_a = g_a + sig_in*(M-0.2*A_neg)*g_max;
    p_a = p_a + A_pos*sig_in;
    
    g_a = g_a + (p_a-0.2*A_pos)*g_max*Po_ex;
    
    g_a = max(g_a,0);
    g_a = min(g_a,g_max);
    
% 	if (Po_ex==1)
% 		M = M-A_neg;
% 	end

% 	count = 0;
% 	count_1 = 0;
% 	for i=1:n
% 		if (sig_in(i)==1)
% % 			count_1 = count_1+1;
% 			g_ex = g_ex + g_a(i);
% 			g_a(i) = g_a(i)+M*g_max; 
% 			p_a(i) = p_a(i)+A_pos;
% 		end
% 
% 		if (Po_ex==1)
% 			g_a(i) = g_a(i)+p_a(i)*g_max;
% 		end
% 
% 		if g_a(i)<0
% 			g_a(i) = 0;
% % 			count = count+1;
% 		elseif g_a(i)>g_max
% 			g_a(i) = g_max;
% 		end
% 	end
    
end
