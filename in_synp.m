function [g_in] = in_synp(g_in_o, sig_in, m, fs);
%g_in_o是上一个时间点的值，sig_in是所有的in-synp的输入信号，m是in-synp的个数
%fs是抽样频率（也就是模拟的最小时间的倒数）

	tao_in = 5.0/1000;
	g_in_ba = 0.05;

	g_in = (-g_in_o)/tao_in*1.0/fs+g_in_o;

	for i=1:m
		if (sig_in(i)==1)
			g_in = g_in + g_in_ba;
		end
	end
end
