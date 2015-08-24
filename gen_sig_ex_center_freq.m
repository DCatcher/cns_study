function [sig_ex_all, sig_in_all] = gen_sig_ex_center_freq(n,m,tmax,fs,lamda,center_freq,mini_freq,maxi_freq)
	time_all = fs*tmax;

	sig_ex_all  = zeros(time_all, n);
	sig_in_all  = zeros(time_all, m);
    
    test_speed = 0;
    
    lamda_new   = lamda*fs;
    
    if length(lamda_new)==1
        lamda_new  = repmat(lamda_new, m, 1);
    else
        disp('error, lambda new should be one number!');
        pause();
    end
    
    pre_all     = cell(n+m, 1);
    
    use_pre     = random('exp', fs*1.0/center_freq, 1, ceil(2*tmax*maxi_freq));
    freq_list   = mini_freq:((maxi_freq - mini_freq)/n):maxi_freq;
    for i=1:n
        pre_all{i}      = use_pre/center_freq*freq_list(i);
        pre_all{i}      = cumsum(pre_all{i});
    end
    
    for i=1:m
        pre_all{n + i}  = random('exp', lamda_new(i), 1, ceil(2*tmax*fs/lamda_new(i)));
        pre_all{n + i}  = cumsum(pre_all{n + i});
    end
    
    tn = tmax*fs;

	for i=1:n
        pre = pre_all{i};
        while (pre(end)<tn)
            disp('error');
            pause();
        end
        ans_lists = find((pre>tn),1);
%         [tmp ans_lists] = findpeaks((pre<tn)*2+(1:length(pre)));
        ans_my = round(pre(1:(ans_lists(1)-1)))+1;        
        sig_ex_all(ans_my,i) = 1;
	end
	for i=1:m
        pre = pre_all{i+n};
        while (pre(end)<tn)
            disp('error');
            pause();
        end        
        ans_lists = find((pre>tn),1);
%         [tmp ans_lists] = findpeaks((pre<tn)*2+(1:length(pre)));        
        ans_my = round(pre(1:ans_lists(1)))+1;          
        sig_in_all(ans_my,i) = 1;
	end
end