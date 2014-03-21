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

function [g_in, M_in, p_a_in, g_a_in] = in_synp(g_in_o, M_in_o, p_a_in_o, g_a_in_o, Po_ex, sig_in, m, fs)

	tao_ex_in = 5.0/1000;
	g_max_in = 0.07;
	tao_neg_in = 30.0/1000;
	tao_pos_in = 30.0/1000;
	A_pos_in = 0.015;
	A_neg_in = A_pos_in;

	g_in = g_in_o+(-g_in_o)/tao_ex_in*1.0/fs;
	M_in = M_in_o + (-M_in_o)/tao_neg_in*1.0/fs;
	p_a_in = p_a_in_o+(-p_a_in_o)/tao_pos_in*1.0/fs;
	g_a_in = g_a_in_o;
	
    M_in = M_in + Po_ex*A_neg_in;
    g_in = g_in + sum(g_a_in.*sig_in);
    g_a_in = g_a_in + sig_in*(M_in-0.4*A_neg_in)*g_max_in;
    p_a_in = p_a_in + A_pos_in*sig_in;
    
    g_a_in = g_a_in + (p_a_in-0.4*A_pos_in)*g_max_in*Po_ex;
    
    g_a_in = max(g_a_in,0);
    g_a_in = min(g_a_in,g_max_in);
    
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
