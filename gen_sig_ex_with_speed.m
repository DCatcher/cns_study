function [sig_ex_all, sig_in_all, speed_ex] = gen_sig_ex_with_speed(n,m,tmax,fs,lamda)
	time_all = fs*tmax;

	sig_ex_all  = zeros(time_all, n);
	sig_in_all  = zeros(time_all, m);
    speed_ex    = zeros(time_all, n);
    
    lamda_new = lamda*fs;
    
    if length(lamda_new)==1
        lamda_new  = repmat(lamda_new, n+m, 1);
    end
    
    pre_all     = cell(n+m, 1);
    
    tn = tmax*fs;
    
    for i=1:(n+m)
        pre_all{i}  = random('exp', lamda_new(i), 1, ceil(2*tmax*fs/lamda_new(i)));
        while sum(pre_all{i})<tn
            pre_all{i}  = [pre_all{i}, random('exp', lamda_new(i), 1, ceil(tmax*fs/lamda_new(i)))];
        end
        pre_all{i}  = cumsum(pre_all{i});
    end

	for i=1:n
        pre = pre_all{i};
%         disp(size(pre));
%         disp(size(pre_all{i}));
        while (pre(end)<tn)
            disp('error');
            disp(pre(end));
            pause();
        end
        ans_lists = find((pre>tn),1);
%         [tmp ans_lists] = findpeaks((pre<tn)*2+(1:length(pre)));
        ans_my = round(pre(1:(ans_lists(1)-1)))+1;        
        sig_ex_all(ans_my,i) = 1;
        sta_time    = 0;
        for j=1:length(ans_my)
            speed_tmp       = 1.2*(100/(ans_my(j) - sta_time));
            speed_ex((sta_time+1):ans_my(j), i)    = speed_tmp;
            
            sta_time        = ans_my(j);
        end
	end
	for i=1:m
        pre = pre_all{i+n};
%         disp(size(pre));
%         disp(size(pre_all{i+n}));
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
