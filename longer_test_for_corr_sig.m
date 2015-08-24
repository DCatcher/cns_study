function ret_val = longer_test_for_corr_sig(cent_freq, anot_freq)

test_num    = 300;

ans_list    = [];

for indx_i=1:test_num
    if mod(indx_i, 40)==0
        fprintf('now indx: %i\n', indx_i);
    end
    ans_list(end+1)     = tmp_test_for_corr_sig(cent_freq, anot_freq, 0);
end

% hist(ans_list);
% disp(mean(ans_list));

ret_val     = mean(ans_list);