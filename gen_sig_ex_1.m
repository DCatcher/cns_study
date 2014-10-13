function [sig_ex_all, sig_in_all] = gen_sig_ex_1(n,m,tmax,fs,lamda)
	time_all = fs*tmax;

	sig_ex_all  = zeros(time_all, n);
	sig_in_all  = zeros(time_all, m);
    
    test_speed = 0;
    
    lamda_new   = lamda*fs;
    
    if length(lamda_new)==1
        lamda_new  = repmat(lamda_new, n+m, 1);
    end
    
    pre_all     = cell(n+m, 1);
    
    for i=1:(n+m)
        pre_all{i}  = random('exp', lamda_new(i), 1, ceil(2*tmax*fs/lamda_new(i)));
        pre_all{i}  = cumsum(pre_all{i});
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