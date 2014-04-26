function [sig_ex_all pic_batch] = gen_sig_pic_temporal(sig_dim,batch_size,time_per,fs,images_in, IMAGES, divide_len, pattern_patches)
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
		sig_pic(:,i) = ceil(tmp(:)*divide_len);
		sig_pic(:,i) = max(sig_pic(:,i), 1);
		sig_pic(:,i) = min(sig_pic(:,i), divide_len);
	end

	for i=1:batch_size
		time_sta = (i-1)*time_per+1;
		time_end = i*time_per;
		for j=1:sig_dim*sig_dim
			sig_ex_all(time_sta:time_end, j) = pattern_patches(sig_pic(j,i), :);
		end
	end
end
