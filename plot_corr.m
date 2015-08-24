function plot_corr(data_path, neuron_num_1, neuron_num_2)

load(data_path);

watch_num   = -60:1:60;
sig_len     = size(sigs, 1);

sig_1   = sigs(:, neuron_num_1);
sig_2   = sigs(:, neuron_num_2);

corr_list   = [];

for indx_i=watch_num
    sig_tmp     = zeros(sig_len, 1);
    if indx_i<0
        sig_tmp(1:(sig_len+indx_i))     = sig_2((1-indx_i):sig_len);
    else
        sig_tmp((1+indx_i):sig_len)     = sig_2(1:(sig_len-indx_i));
    end
    
    corr_list(end+1)    = sum(sig_1.*sig_tmp);
end

disp(corr_list);

plot(watch_num, corr_list);