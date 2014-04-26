function [sig_ex_all,pic_batch] = gen_sig_pic(sig_dim,batch_size,time_per,fs,images_in, IMAGES, lamda_max, lamda_min)
	time_all = time_per*batch_size;

	sig_ex_all = zeros(time_all, sig_dim*sig_dim);
	sig_pic = zeros(sig_dim*sig_dim,batch_size);

	[x y z] = size(images_in);
	which_x = ceil(rand(batch_size,1)*(x-sig_dim));
	which_y = ceil(rand(batch_size,1)*(y-sig_dim));
	which_z = ceil(rand(batch_size,1)*z);
	pic_batch = zeros(batch_size, sig_dim, sig_dim);

	for i =1:batch_size
% 		tmp = images_in((which_x(i)):(which_x(i)+sig_dim-1),(which_y(i)):(which_y(i)+sig_dim-1),which_z(i));
		pic_batch(i, :, :) = IMAGES((which_x(i)):(which_x(i)+sig_dim-1),(which_y(i)):(which_y(i)+sig_dim-1),which_z(i));
		tmp = pic_batch(i,:,:);
		tmp = tmp-mean2(tmp);
		tmp = tmp/max(max(max(abs(tmp))));
		tmp = tmp/2+0.5;
		sig_pic(:,i) = tmp(:)*(lamda_max-lamda_min)+lamda_min;
	end
	r_a = sig_pic;

	time_a = time_per.*r_a/fs;
	max_time = max(max(time_a));
	max_time = ceil(max_time)+1;
	max_time = max_time*4;
	
	lamda_standrad = 0.1*fs;
	big_exp_rand = random('exp', lamda_standrad, sig_dim*sig_dim, batch_size, max_time);
	r_a_new = zeros(sig_dim*sig_dim,batch_size,max_time);
	for i=1:max_time
		r_a_new(:,:,i) = r_a;
	end
	big_exp_rand = big_exp_rand.*(fs./r_a_new)/lamda_standrad;
	interval_my = ones(1, batch_size)*time_per;

	for i=1:(sig_dim*sig_dim)
		middle_rand = zeros(batch_size, max_time);
		middle_rand(:,:) = big_exp_rand(i,:,:);
		middle_rand = cumsum(middle_rand')';
		new_start = [1 interval_my(1:(end-1))];
		new_start = cumsum(new_start);
		middle_rand = middle_rand + repmat(new_start', 1, max_time);
		new_end = cumsum(interval_my)+1;
		ans_my = find(middle_rand<repmat(new_end',1, max_time));
		sig_ex_all(round(middle_rand(ans_my)), i) = 1;
	end
end
