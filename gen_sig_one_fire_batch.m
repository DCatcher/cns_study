function [sig_ex_all, pic_batch]=gen_sig_one_fire_batch(sig_dim, batch_size, time_per, fs, IMAGES, fire_interval)
	time_all = time_per*batch_size+1;
	sta_time = 1;

	sig_ex_all = zeros(time_all, sig_dim^2);

	[x,y,z] = size(IMAGES);
	which_x = ceil(rand(batch_size,1)*(x-sig_dim));
	which_y = ceil(rand(batch_size,1)*(y-sig_dim));
	which_z = ceil(rand(batch_size,1)*z);
	pic_batch = zeros(batch_size, sig_dim, sig_dim);

	for i =1:batch_size
% 		tmp = images_in((which_x(i)):(which_x(i)+sig_dim-1),(which_y(i)):(which_y(i)+sig_dim-1),which_z(i));
		pic_batch(i, :, :) = IMAGES((which_x(i)):(which_x(i)+sig_dim-1),(which_y(i)):(which_y(i)+sig_dim-1),which_z(i));
		pattern_fre = pic_batch(i,:,:);

		[tmp,indx] = sort(pattern_fre(:),'descend');

        add_time = 1;
        inc_time = 0;
		for j=1:sig_dim^2
            if tmp(j)>0
                sig_ex_all(sta_time+round(add_time),indx(j)) = 1;
                add_time = add_time+2+inc_time;
                inc_time = inc_time + 0.5;
            end
		end
		sta_time = sta_time+time_per;
	end
end
