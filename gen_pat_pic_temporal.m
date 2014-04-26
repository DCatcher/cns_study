function [pattern_patches]=gen_pat_pic_temporal(param)
	fre_min 		= param.fre_min;
	fre_max 		= param.fre_max; 
	time_per 		= param.time_per; 
	divide_len 		= param.divide_len;
	jitter_stride 	= param.jitter_stride;
	fs				= param.fs;

	pattern_patches = zeros(divide_len, time_per);
	pre_exp = random('exp', fs/fre_min, 1, round(4*fre_min*time_per/fs));
	sum_exp = cumsum(pre_exp);
	ans_exp = sum_exp(find(sum_exp<time_per));
	pattern_patches(1,round(ans_exp)+1) = 1;

	fre_stride = (fre_max - fre_min)/divide_len/fs*time_per;

	for i=2:divide_len

		while sum(pattern_patches(i,:))< sum(pattern_patches(i-1,:))
			pattern_patches(i,:) = zeros(1, time_per);
			tmp_exp = find(pattern_patches(i-1,:));
			tmp_exp = tmp_exp + randn(size(tmp_exp))*jitter_stride;

			tmp_exp = min(tmp_exp, time_per);
			tmp_exp = max(tmp_exp, 1);
			pattern_patches(i,round(tmp_exp)) = 1;
		end

		if rand < fre_stride
			tmp_idx = ceil(rand*time_per);
			while pattern_patches(i, tmp_idx)==1
				tmp_idx = ceil(rand*time_per);
			end
			pattern_patches(i, tmp_idx)=1;
		end
	end
end
