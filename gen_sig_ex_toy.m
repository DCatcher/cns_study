function [sig_ex_all, sig_in_all] = gen_sig_ex_toy(n,m,tmax,fs,lamda_ex, lamda_in);
	time_all = fs*tmax;

	sig_ex_all = zeros(time_all, n);
	sig_in_all = zeros(time_all, m);
    
    lamda_new = fs;
    lamda = min(min(lamda_ex), min(lamda_in));
    
    pre_all = random('exp', lamda_new, n+m, ceil(4*tmax/lamda));
    pre_all = cumsum(pre_all')';
    
    tn = tmax*fs;
    t_left = 1+n+m;

	for i=1:n
        pre = pre_all(i,:)*lamda_ex(i);
        while (pre(end)<tn)
            pre = [pre pre_all(t_left,:)+pre(end)];
            t_left = t_left+1;
        end
        ans_lists = find((pre>tn),1);
        ans_my = round(pre(1:(ans_lists(1)-1)))+1;        
        sig_ex_all(ans_my,i) = 1;
	end
	for i=1:m
        pre = pre_all(i+n,:)*lamda_in(i);
        while (pre(end)<tn)
            pre = [pre pre_all(t_left,:)+pre(end)];
            t_left = t_left+1;
        end        
        ans_lists = find((pre>tn),1);      
        ans_my = round(pre(1:ans_lists(1)))+1;          
        sig_in_all(ans_my,i) = 1;
	end
end
