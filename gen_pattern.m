function pattern = gen_pattern(time_pa, fs, wid_pa, lamda)
	time_len = time_pa*fs;
	lamda_new = lamda*fs;

	pre_all = random('exp', lamda_new, wid_pa+time_pa/lamda, time_pa/lamda+20);
	pre_all = cumsum(pre_all')';
	t_left = 1+wid_pa;
	pattern = zeros(time_len, wid_pa);

	for i=1:wid_pa
        pre = pre_all(i,:);
        while (pre(end)<time_len)
            pre = [pre pre_all(t_left,:)+pre(end)];
            t_left = t_left+1;
        end
        ans_lists = find((pre>time_len),1);
        ans_my = round(pre(1:(ans_lists(1)-1)))+1;        
        pattern(ans_my,i) = 1;
	end
end
