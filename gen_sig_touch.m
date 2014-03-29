function [sig_ex_all, sig_in_all] = gen_sig_touch(n,m,tmax,fs,lamda,r_1,tao_c,sigma);
	time_all = fs*tmax;

	sig_ex_all = zeros(time_all, n);
	sig_in_all = zeros(time_all, m);
    
    test_speed = 0;
    
    lamda_new = lamda*fs;
	tao_c_new = tao_c*fs;
    
    pre_all = random('exp', lamda_new, m+tmax/lamda, round(tmax/lamda*2));
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

    pre_all = random('exp', tao_c_new, 2, round(tmax/tao_c*2));
    pre_sum = cumsum(pre_all')';
    
    t_left = 2;

	pre = pre_sum(1,:);
	interval_my = pre_all(1,:);
	while (pre(end)<tn)
		pre = [pre pre_sum(t_left,:)+pre(end)];
		interval_my = [interval_my pre_all(t_left,:)];
		t_left = t_left+1;
	end        
	ans_lists = find((pre>tn),1);
	interval_my = interval_my(1:ans_lists(1));

	len_inter = length(interval_my);
	s = rand(1, len_inter)*n;
	s = repmat(s,n,1);

	r_ba = 1/lamda;
	a = (1:n)';
	a = repmat(a,1,len_inter);
	r_a = r_ba + r_1*(exp(-(s-a).^2/(2*sigma^2))+exp(-(s-a+n).^2/(2*sigma^2))+exp(-(s-a-n).^2/(2*sigma^2)));
	r_a = max(r_a, 0.001);

	time_a = repmat(interval_my, n, 1).*r_a/fs;
	max_time = max(max(time_a));
	max_time = ceil(max_time)+1;
	max_time = max_time*2;
	
	lamda_standrad = 0.1*fs;
	big_exp_rand = random('exp', lamda_standrad, n, len_inter, 2*max_time);
	big_exp_rand = big_exp_rand.*(fs./repmat(r_a,1,1,2*max_time))/lamda_standrad;

	for i=1:n
		middle_rand = zeros(len_inter, 2*max_time);
		middle_rand(:,:) = big_exp_rand(i,:,:);
		middle_rand = cumsum(middle_rand')';
		new_start = [1 interval_my(1:(end-1))];
		new_start = cumsum(new_start);
		middle_rand = middle_rand + repmat(new_start', 1, 2*max_time);
		new_end = cumsum(interval_my)+1;
		ans_my = find(middle_rand<repmat(new_end',1,2*max_time));
		sig_ex_all(round(middle_rand(ans_my)), i) = 1;
	end
end
