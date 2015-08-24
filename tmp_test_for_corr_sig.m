function ret_val    = tmp_test_for_corr_sig(cent_freq, anot_freq, show_flag, test_flag)

if nargin < 1
    cent_freq   = 50;
    anot_freq   = 51;
end

if nargin < 3
    show_flag   = 1;
end

if nargin < 4
    test_flag   = 0;
end

pre_1   = random('exp', 1.0/cent_freq * 1000, 1, 1000);
% pre_2   = pre_1/cent_freq*anot_freq;
pre_2   = random('exp', 1.0/anot_freq * 1000, 1, 1000);

sig_pos_1   = cumsum(pre_1);
sig_pos_2   = cumsum(pre_2);

time_len    = 10000;

sig_1   = zeros(1, time_len);
sig_2   = zeros(1, time_len);

sig_pos_1   = sig_pos_1(sig_pos_1 < length(sig_1));
sig_pos_2   = sig_pos_2(sig_pos_2 < length(sig_2));
% disp(sig_pos_1);

sig_1(floor(sig_pos_1) + 1)         = 1;
sig_2(floor(sig_pos_2) + 1)         = 1;

if test_flag==1
    sig_1   = (sig_1 - mean(sig_1))/std(sig_1);
    sig_2   = (sig_2 - mean(sig_2))/std(sig_2);
end

watch_num   = -60:1:60;
sig_len     = size(sig_1, 2);
% 
% sig_1   = sigs(:, neuron_num_1);
% sig_2   = sigs(:, neuron_num_2);

corr_list   = [];

for indx_i=watch_num
    sig_tmp     = zeros(1, sig_len);
    if indx_i<0
        sig_tmp(1:(sig_len+indx_i))     = sig_2((1-indx_i):sig_len);
    else
        sig_tmp((1+indx_i):sig_len)     = sig_2(1:(sig_len-indx_i));
    end
    
    corr_list(end+1)    = sum(sig_1.*sig_tmp);
end

if test_flag==0

    ret_val     = sum(corr_list.*sign(watch_num));

    if show_flag==1
        plot(watch_num, corr_list);
        disp(ret_val);
    end
else
    ret_val     = corr_list;
    if show_flag==1
        plot(watch_num, corr_list);
    end
end