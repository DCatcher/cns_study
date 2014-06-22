
fs = 1000;
tao_neg_in = 0.02;
tao_pos_in = 0.02;
g_max_in = 0.015;
A_neg_in = 0.005;
A_pos_in = 0.005;

time_vec = -50:50;
delta_g = [];

for dtime=time_vec
	if dtime<0
		sig_in = [1 zeros(1,abs(dtime))];
		Po_ex = [zeros(1,abs(dtime)) 1];
	else
		Po_ex = [1 zeros(1,abs(dtime))];
		sig_in = [zeros(1,abs(dtime)) 1];
	end

	g_a_in = 0.5*g_max_in;
	M_in = 0;
	p_a_in = 0;
	for i=1:(abs(dtime)+1)
		M_in = M_in + (-M_in)/tao_neg_in*1.0/fs;
		p_a_in = p_a_in+(-p_a_in)/tao_pos_in*1.0/fs;

		M_in = M_in + Po_ex(i)*A_neg_in;
		g_a_in = g_a_in + sig_in(i)*(M_in-0.4*A_neg_in)*g_max_in;
		p_a_in = p_a_in + A_pos_in*sig_in(i);

		g_a_in = g_a_in + (p_a_in-0.4*A_pos_in)*g_max_in*Po_ex(i);

		g_a_in = max(g_a_in,0);
		g_a_in = min(g_a_in,g_max_in);
	end

	delta_g(end+1) = g_a_in/g_max_in-0.5;
end

plot(time_vec, delta_g);
