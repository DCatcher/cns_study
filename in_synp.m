function [g_in] = in_synp(g_in_o, sig_in, m, fs);
%g_in_o����һ��ʱ����ֵ��sig_in�����е�in-synp�������źţ�m��in-synp�ĸ���
%fs�ǳ���Ƶ�ʣ�Ҳ����ģ�����Сʱ��ĵ�����

	tao_in = 5.0/1000;
	g_in_ba = 0.05;

	g_in = (-g_in_o)/tao_in*1.0/fs+g_in_o;

	for i=1:m
		if (sig_in(i)==1)
			g_in = g_in + g_in_ba;
		end
	end
end
