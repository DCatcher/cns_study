function [sig_ex_all,sig_in_all] = gen_sig_one_fire(n,m,tmax,fs,pattern_fre,fire_interval)
	time_all = fs*tmax;

	sig_ex_all = zeros(time_all, n);
	sig_in_all = zeros(time_all, m);

	[tmp indx] = sort(pattern_fre,'descend');

	for i=1:n
		sig_ex_all(find(indx==i)*fire_interval,i) = 1;
	end
    lamda_new = 0.1*fs;
	lamda = 0.1;
    
    pre_all = random('exp', lamda_new, m+tmax/lamda, tmax/lamda+20);
    pre_all = cumsum(pre_all')';
    
    tn = tmax*fs;
    t_left = 1+m;
	for i=1:m
        pre = pre_all(n,:);
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
