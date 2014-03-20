function [sig_ex_all, sig_in_all] = gen_sig_ex_1(n,m,tmax,fs,lamda);
	time_all = fs*tmax;

	sig_ex_all = zeros(time_all, n);
	sig_in_all = zeros(time_all, m);
    
    test_speed = 0;
    
    lamda_new = lamda*fs;
    
    pre_all = random('exp', lamda_new, n+m+tmax/lamda, tmax/lamda+20);
    pre_all = cumsum(pre_all')';
    
    tn = tmax*fs;
    t_left = 1+n+m;

	for i=1:n
        pre = pre_all(i,:);
        while (pre(end)<tn)
            pre = [pre pre_all(t_left,:)+pre(end)];
            t_left = t_left+1;
        end
        ans_lists = find((pre>tn),1);
%         [tmp ans_lists] = findpeaks((pre<tn)*2+(1:length(pre)));
        ans_my = round(pre(1:ans_lists(1)))+1;        
        sig_ex_all(ans_my,i) = 1;
	end
	for i=1:m
        pre = pre_all(i+n,:);
        while (pre(end)<tn)
            pre = [pre pre_all(t_left,:)+pre(end)];
            t_left = t_left+1;
        end        
        ans_lists = find((pre>tn),1);
%         [tmp ans_lists] = findpeaks((pre<tn)*2+(1:length(pre)));        
        ans_my = round(pre(1:ans_lists(1)))+1;          
        sig_in_all(ans_my,i) = 1;
	end
end
