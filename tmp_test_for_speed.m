%% For speed test

tmax        = 10;
freq        = 1000;

fire_freq   = 40;
fire_list   = double(rand(tmax*freq, 1) < fire_freq/freq);

sig_freq        = 20;
sig_in_list     = double(rand(tmax*freq, 1) < sig_freq/freq);

stdp_sig_to_ex_neg  = 0:60;
num_sig_to_ex_neg   = zeros(length(stdp_sig_to_ex_neg), 1);

stdp_sig_to_ex_pos  = -60:-1;
num_sig_to_ex_pos   = zeros(length(stdp_sig_to_ex_pos), 1);

rep_time        = 1000;

tic
for indx_i=1:rep_time
    tmp_fire_list       = find(fire_list(:)==1);
    tmp_sig_in_list     = find(sig_in_list(:)==1);

    for k=1:length(tmp_fire_list)
        tmp_fire_time                       = tmp_fire_list(k);
        tmp_sig_time_list                   = (tmp_sig_in_list >= tmp_fire_time) & (tmp_sig_in_list < (tmp_fire_time + length(stdp_sig_to_ex_neg) - 1));
    %     param.delta_g_a_sig_to_ex(i,j)      = param.delta_g_a_sig_to_ex(i,j) + sum(param.stdp_sig_to_ex_neg(tmp_sig_in_list(tmp_sig_time_list) - tmp_fire_time + 1));
    %     param.delta_t_his_sig_to_ex_neg{i,j}(tmp_sig_in_list(tmp_sig_time_list) - tmp_fire_time + 1) ...
    %                                         = param.delta_t_his_sig_to_ex_neg{i,j}(tmp_sig_in_list(tmp_sig_time_list) - tmp_fire_time + 1) + 1;
        num_sig_to_ex_neg(tmp_sig_in_list(tmp_sig_time_list) - tmp_fire_time + 1)   = num_sig_to_ex_neg(tmp_sig_in_list(tmp_sig_time_list) - tmp_fire_time + 1) + 1;
        tmp_sig_time_list                   = (tmp_sig_in_list < tmp_fire_time) & (tmp_sig_in_list > (tmp_fire_time - length(stdp_sig_to_ex_pos)));
    %     param.delta_g_a_sig_to_ex(i,j)      = param.delta_g_a_sig_to_ex(i,j) + sum(param.stdp_sig_to_ex_pos(tmp_fire_time - tmp_sig_in_list(tmp_sig_time_list)));
    %     param.delta_t_his_sig_to_ex_pos{i,j}(tmp_fire_time - tmp_sig_in_list(tmp_sig_time_list)) ...
    %                                         = param.delta_t_his_sig_to_ex_pos{i,j}(tmp_fire_time - tmp_sig_in_list(tmp_sig_time_list)) + 1;
        num_sig_to_ex_pos(tmp_fire_time - tmp_sig_in_list(tmp_sig_time_list))       = num_sig_to_ex_pos(tmp_fire_time - tmp_sig_in_list(tmp_sig_time_list)) + 1;
    end
end
toc

tic
for indx_i=1:rep_time
    corr_list   = [];
    for indx_j=stdp_sig_to_ex_neg
        corr_list(end+1)    = sum(fire_list(1:(length(fire_list) - indx_j)).*sig_in_list((1+indx_j):length(fire_list)));
    end
    for indx_j=stdp_sig_to_ex_pos
        corr_list(end+1)    = sum(fire_list((1 - indx_j):length(fire_list)).*sig_in_list(1:(length(fire_list) + indx_j)));
    end
end
toc
