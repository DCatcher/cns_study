function [sig_ex_all, sig_in_all] = gen_sig_ex_corr(n,m,tmax,fs,lamda,r_a);
	time_all = fs*tmax;

	sig_ex_all = zeros(time_all, n);
	sig_in_all = zeros(time_all, m);
    
    test_speed = 0;
    
    lamda_new = lamda*fs;
    
    pre_all = random('exp', lamda_new, m+tmax/lamda, tmax/lamda+20);
    pre_all = cumsum(pre_all')';
    
    tn = tmax*fs;
    t_left = 1+m;

	for i=1:m
        pre = pre_all(i,:);
        while (pre(end)<tn)
            pre = [pre pre_all(t_left,:)+pre(end)];
            t_left = t_left+1;
        end        
        ans_lists = find((pre>tn),1);
        ans_my = round(pre(1:ans_lists(1)))+1;          
        sig_in_all(ans_my,i) = 1;
	end

	r_max = max(r_a);
	pre_all = rand(n, round(r_max*tmax*2));
	tmp_mat = repmat(r_a', 1, round(r_max*tmax*2));
	pre_all = -log(1-pre_all).*(fs./tmp_mat);
	pre_all = cumsum(pre_all')';
    t_left = 1+n;
	
	for i=1:n
        pre = pre_all(i,:);
        while (pre(end)<tn)
            pre = [pre pre_all(t_left,:)+pre(end)];
            t_left = t_left+1;
        end        
        ans_lists = find((pre>tn),1);
        ans_my = round(pre(1:ans_lists(1)))+1;          
        sig_ex_all(ans_my,i) = 1;
	end
end
