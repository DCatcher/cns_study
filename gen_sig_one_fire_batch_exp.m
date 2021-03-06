function [sig_ex_all, pic_batch]=gen_sig_one_fire_batch_exp(sig_dim, batch_size, time_per, IMAGES, fire_max_time, I_0, param)
	time_all = time_per*batch_size+1;
	sta_time = 1;

	sig_ex_all = zeros(time_all, 2*sig_dim^2);

	[x,y,z] = size(IMAGES);
	which_x = ceil(rand(batch_size,1)*(x-sig_dim));
	which_y = ceil(rand(batch_size,1)*(y-sig_dim));
	which_z = ceil(rand(batch_size,1)*z);
	pic_batch = zeros(batch_size, sig_dim, sig_dim);
    add_time_list = [];
	for i =1:batch_size
% 		tmp = images_in((which_x(i)):(which_x(i)+sig_dim-1),(which_y(i)):(which_y(i)+sig_dim-1),which_z(i));
		pic_batch(i, :, :) = IMAGES((which_x(i)):(which_x(i)+sig_dim-1),(which_y(i)):(which_y(i)+sig_dim-1),which_z(i));
		pattern_fre = pic_batch(i,:,:);
        
		for j=1:sig_dim^2
            if pattern_fre(j)>0
%                 add_time = fire_max_time*exp(-pattern_fre(j)/I_0);
                add_time = param.fire_base_time_ex + fire_max_time*(1-(1/(1+(exp(-pattern_fre(j)/I_0)))-0.5)*2);
                add_time_list(end+1) = add_time;
%                 add_time = add_time:add_time:fire_max_time;
                tmp_add  = 3;
                if (pattern_fre(j)>2*I_0) 
                    tmp_add = tmp_add +1;
                end
                
                if (pattern_fre(j)>3*I_0) 
                    tmp_add = tmp_add +2;
                end
                
                if (pattern_fre(j)>4*I_0) 
                    tmp_add = tmp_add +2;
                end
                
                add_time = add_time:1:min(add_time+tmp_add, time_per);
                sig_ex_all(sta_time+round(add_time),j) = 1;
            else
%                 add_time = fire_max_time*exp(pattern_fre(j)/I_0);
                add_time = param.fire_base_time_in + fire_max_time*(1-(1/(1+(exp( pattern_fre(j)/I_0)))-0.5)*2);
                add_time_list(end+1) = add_time;
%                 add_time = add_time:add_time:fire_max_time;
%                 add_time = add_time:1:min(add_time+1, time_per);

                tmp_add  = 3;
                if (abs(pattern_fre(j))>2*I_0) 
                    tmp_add = tmp_add +1;
                end
                
                if (abs(pattern_fre(j))>3*I_0) 
                    tmp_add = tmp_add +2;
                end
                
                if (abs(pattern_fre(j))>4*I_0) 
                    tmp_add = tmp_add +2;
                end
                
                add_time = add_time:1:min(add_time+tmp_add, time_per);
                sig_ex_all(sta_time+round(add_time),j+sig_dim^2) = 1;
            end
        end
        
		sta_time = sta_time+time_per;
    end
%     hist(add_time_list,0:1:100);
%     disp(add_time_list(1:(64*1)));
    pause(1);
end
