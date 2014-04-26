function neuron_sparse_coding(param)
	sig_dim = param.sig_dim;
	ex_n = param.ex_n;
	fs = param.fs;
	train_num = param.train_num;
	lamda_max = param.lamda_max;
    lamda_min = param.lamda_min;
	time_per = param.time_per;
	batch_size = param.batch_size;
	tr_mat = param.tr_mat;

    cont = param.cont;
    
	tao_m = param.tao_m;
	v_rest = param.v_rest;
	e_ex = param.e_ex;
	e_in = param.e_in;
	v_th = param.v_th;
	v_reset = param.v_reset;
    
	tao_ex = param.tao_ex;
	g_max = param.g_max;
	tao_neg = param.tao_neg;
	tao_pos = param.tao_pos;
	A_pos = param.A_pos;
	A_neg = param.A_neg;
    
    tao_ex_in = param.tao_ex_in;
	g_max_in = param.g_max_in;
	tao_neg_in = param.tao_neg_in;
	tao_pos_in = param.tao_pos_in;
	A_pos_in = param.A_pos_in;
	A_neg_in = param.A_neg_in;
	neg_amp = param.neg_amp;
    
    short_report_mode = param.short_report_mode;
	display_mode = param.display_mode;
    batch_learn = param.batch_learn;
    easy_inhib = param.easy_inhib;
	time_pattern = param.time_pattern;
	divide_len = param.divide_len;

	load(tr_mat);
	images_in = IMAGES-mean2(IMAGES); 
	images_in = images_in/std2(images_in);
%	images_in = 1./(1+exp(images_in));
	images_in = images_in/max(max(max(images_in)))/2+0.5;
	images_in = images_in*lamda_max;

    time_all = time_per*batch_size;
	sig_n = sig_dim*sig_dim;
    vol_ex = zeros(ex_n,time_all);
	vol_ex(:,1) = v_reset;
	sta_pic = zeros(ex_n, sig_dim, sig_dim);
    sta_pic_ans = zeros(ex_n, sig_dim, sig_dim);
	sta_num = zeros(ex_n, 1);

	sta_figure = figure;
	spike_figure = figure;
	sper_patch_figure = figure;
	for lamda_now = lamda_max
        tag_now = ex_n;
		if (cont==0)
			g_a_sig_to_ex_use = g_max*rand(ex_n, sig_n);
			p_a_sig_to_ex = zeros(ex_n, sig_n);
			M_sig_to_ex = 0*ones(ex_n,1);
			g_sig_to_ex = 1*ones(ex_n,1);
			
			g_a_ex_to_other_use = g_max_in*rand(ex_n, ex_n);
			p_a_ex_to_other = zeros(ex_n, ex_n);
			M_ex_to_other = 0*ones(ex_n,1);
			g_ex_to_other = 0*ones(ex_n,1);
		elseif (cont==1)
			load save_data
        end

        batch_num_count = 0;
		if time_pattern==1
			pattern_patches = gen_pat_pic_temporal(param);
		end
		for patch_now=0:batch_size:train_num
			spike_per_patch = zeros(1, batch_size);

            if short_report_mode==1
                fprintf('patch_now:%i/%i, step: %i\n',patch_now,train_num,batch_size);
            end
			if time_pattern==0
				[sig_ex_all,pic_batch] = gen_sig_pic(sig_dim,batch_size,time_per,fs,images_in,IMAGES, lamda_max, lamda_min);
			else
				[sig_ex_all,pic_batch] = gen_sig_pic_temporal(sig_dim,batch_size,time_per,fs,images_in, IMAGES, divide_len, pattern_patches);
			end
			if (patch_now>0)
				vol_ex(:,1) = vol_ex(:,time_all);
			end
			batch_indx = 1;
			batch_fire_num = zeros(ex_n, batch_size);
            g_a_sig_to_ex = g_a_sig_to_ex_use;
            g_a_ex_to_other = g_a_ex_to_other_use;
			for i=2:time_all
				if mod(i, time_per)==1
					batch_fire_num(:,batch_indx) = sum((vol_ex(:,(i-time_per):(i-1))>-10)')';
					spike_per_patch(batch_indx) = sum(batch_fire_num(:,batch_indx)>0);
					for j=1:ex_n
						sta_pic(j,:,:) = batch_fire_num(j,batch_indx)*pic_batch(batch_indx,:,:)+sta_pic(j,:,:);
						sta_num(j) = sta_num(j)+batch_fire_num(j,batch_indx);
					end
					batch_indx = batch_indx+1;
					vol_ex(:,i-1) = v_reset;
				end
                if short_report_mode==1
                    step_size = 100;
                    if (mod(i*step_size,time_all)==(-1))
                        fprintf('subpatch_now:%i/%i, excited:%i, sig:%i\n',floor(i*step_size/time_all),step_size,sum(vol_ex(:,i-1)>-10), sum(sig_ex_all(i,:)));
                    end
                end                
				vol_ex(:,i) = v_reset*(vol_ex(:,i-1)>v_th);
                vol_ex(:,i) = vol_ex(:,i) + (vol_ex(:,i-1)<v_th).*...
                    ((v_rest-vol_ex(:,i-1)+g_sig_to_ex.*(e_ex-vol_ex(:,i-1))+g_ex_to_other.*(e_in-vol_ex(:,i-1)))/tao_m*1.0/fs+vol_ex(:,i-1));

				Po_ex = (vol_ex(:,i)>=v_th);
				Po_ex_row = find(Po_ex==1);
				vol_ex(:,i) = vol_ex(:,i).*(vol_ex(:,i)<v_th); 
%				size(sta_pic(Po_ex_row,:,:))
%				size(repmat(pic_batch(floor((i-2)/50)+1,:,:), length(Po_ex_row), 1))
				%{
				for j = 1:min(i-2,10)
					sta_pic(Po_ex_row,:,:) = sta_pic(Po_ex_row,:,:)+repmat(pic_batch(floor((i-2-j)/50)+1,:,:), length(Po_ex_row), 1);
					sta_num(Po_ex_row) = sta_num(Po_ex_row)+1;
				end
				%}
%sig_to_ex begin
				sig_in = sig_ex_all(i,:);
				sig_in_row = find(sig_in==1);

				g_sig_to_ex = g_sig_to_ex*(1-1/tao_ex*1.0/fs);
				M_sig_to_ex = M_sig_to_ex*(1-1/tao_neg*1.0/fs);
				p_a_sig_to_ex = p_a_sig_to_ex*(1-1/tao_pos*1.0/fs);

				M_sig_to_ex = M_sig_to_ex - Po_ex*A_neg;
				if batch_learn==1
                    g_sig_to_ex = g_sig_to_ex + (sum((g_a_sig_to_ex_use(:,sig_in_row))'))';
                else
                    g_sig_to_ex = g_sig_to_ex + (sum((g_a_sig_to_ex(:,sig_in_row))'))';
                end
                
				g_a_sig_to_ex(:,sig_in_row) = g_a_sig_to_ex(:,sig_in_row) + repmat(M_sig_to_ex, 1, length(sig_in_row))*g_max;
				g_a_sig_to_ex(:,sig_in_row) = max(g_a_sig_to_ex(:,sig_in_row),0);
				p_a_sig_to_ex(:,sig_in_row) = p_a_sig_to_ex(:,sig_in_row) + A_pos;

				if length(Po_ex_row)<0.5*ex_n
					g_a_sig_to_ex(Po_ex_row,:) = g_a_sig_to_ex(Po_ex_row,:) + p_a_sig_to_ex(Po_ex_row,:)*g_max;
					g_a_sig_to_ex(Po_ex_row,:) = min(g_a_sig_to_ex(Po_ex_row,:),g_max);
				else
					g_a_sig_to_ex = g_a_sig_to_ex+p_a_sig_to_ex*g_max;
					Po_ex_row = find(Po_ex==0);
					g_a_sig_to_ex(Po_ex_row,:) = g_a_sig_to_ex(Po_ex_row,:) - p_a_sig_to_ex(Po_ex_row,:)*g_max;
					g_a_sig_to_ex = min(g_a_sig_to_ex,g_max);
					Po_ex_row = find(Po_ex==1);
				end
%sig_to_ex finish

%ex_to_other begin
				Po_all = Po_ex;

				sig_in = vol_ex(:,i)'>v_th;
				sig_in_row = find(sig_in==1);

				g_ex_to_other = g_ex_to_other*(1-1/(tao_ex_in*fs));
				if easy_inhib==1
					last_fire_num = sum(vol_ex(:,i-1)>-10);
					g_ex_to_other = g_ex_to_other + (last_fire_num - (vol_ex(:,i-1)>-10))*param.inhibi_strength;
				else
					M_ex_to_other = M_ex_to_other*(1-1/(tao_neg_in*fs));
					p_a_ex_to_other = p_a_ex_to_other*(1-1/(tao_pos_in*fs));

					M_ex_to_other = M_ex_to_other + Po_all*A_neg_in;
					if batch_learn==1
						g_ex_to_other = g_ex_to_other + (sum((g_a_ex_to_other_use(:,sig_in_row))'))';
					else
						g_ex_to_other = g_ex_to_other + (sum((g_a_ex_to_other(:,sig_in_row))'))';
					end
					
					g_a_ex_to_other(:,sig_in_row) = g_a_ex_to_other(:,sig_in_row) + repmat(M_ex_to_other-neg_amp*A_neg_in,1,length(sig_in_row))*g_max_in;
					g_a_ex_to_other(:,sig_in_row) = max(g_a_ex_to_other(:,sig_in_row),0);
					g_a_ex_to_other(:,sig_in_row) = min(g_a_ex_to_other(:,sig_in_row),g_max_in);
					p_a_ex_to_other(:,sig_in_row) = p_a_ex_to_other(:,sig_in_row) + A_pos_in;

					if length(Po_ex_row)<0.5*ex_n
						g_a_ex_to_other(Po_ex_row,:) = g_a_ex_to_other(Po_ex_row,:) + (p_a_ex_to_other(Po_ex_row,:)-neg_amp*A_pos_in)*g_max_in;
						g_a_ex_to_other(Po_ex_row,:) = min(g_a_ex_to_other(Po_ex_row,:),g_max_in);
						g_a_ex_to_other(Po_ex_row,:) = max(g_a_ex_to_other(Po_ex_row,:),0);
					else
						g_a_ex_to_other = g_a_ex_to_other+(p_a_ex_to_other-neg_amp*A_pos_in)*g_max_in;
						Po_ex_row = find(Po_ex==0);
						g_a_ex_to_other(Po_ex_row,:) = g_a_ex_to_other(Po_ex_row,:) - (p_a_ex_to_other(Po_ex_row,:)-neg_amp*A_pos_in)*g_max_in;
						g_a_ex_to_other = min(g_a_ex_to_other,g_max_in);
						g_a_ex_to_other = max(g_a_ex_to_other,0);
						Po_ex_row = find(Po_ex==1);
					end
				end
%in_to_other finish
            end
            
            if batch_learn==1
                g_a_sig_to_ex_use = g_a_sig_to_ex;
                g_a_ex_to_other_use = g_a_ex_to_other;
            end
            
			spike_every = sum(([vol_ex']>-10));
			ans_my = mean(sum(([vol_ex']>-10)));
            if short_report_mode==1
                fprintf('spiking_time:%f\n',ans_my);
            end
			if display_mode==1
				figure(sta_figure);
				now_indx = 1;
				sqrt_ex_n = round(sqrt(ex_n));
                big_pic = zeros(sqrt_ex_n*sig_dim, sqrt_ex_n*sig_dim);
                inter_val = 1.0/sqrt_ex_n;
				for i=1:sqrt_ex_n
					for j=1:sqrt_ex_n
                        subplot('Position',[(i-1)*inter_val,(j-1)*inter_val,inter_val*0.9,inter_val*0.9]);
						tmp = zeros(sig_dim, sig_dim);
						tmp(:,:) = sta_pic_ans(now_indx,:,:);
                        tmp = tmp-mean2(tmp);
						tmp = tmp/max(max(abs(tmp))+0.0001)/2+0.5;
						imagesc(tmp,[0 1]);
						axis image off; colormap gray;
						now_indx = now_indx+1;
                        big_pic(((i-1)*sig_dim+1):((i-1)*sig_dim+sig_dim),((j-1)*sig_dim+1):((j-1)*sig_dim+sig_dim)) = tmp;
					end
                end
                %imshow(big_pic/g_max,[0 1]);
				pause(0.1);
				figure(spike_figure);
				hist(spike_every);
				figure(sper_patch_figure);
				hist(spike_per_patch);
            end
%             batch_num_count = batch_num_count+1;
%             if batch_num_count==5
%                 sta_pic = zeros(ex_n, sig_dim, sig_dim);
%                 sta_num = zeros(ex_n, 1);
%                 batch_num_count = 0;
%             end
            
            tmp_sta = find(sta_num>param.sta_big_num);
            for tmp_idx=tmp_sta'
                sta_pic_ans(tmp_idx,:,:) = sta_pic(tmp_idx,:,:)./sta_num(tmp_idx);
                sta_num(tmp_idx) = 0;
                sta_pic(tmp_idx,:,:) = zeros(1, sig_dim, sig_dim);
            end
		end
		sig_ex_all = [];
        s = ['data_',int2str(tag_now)];
        save(s);
        fprintf('saved!%s\n',s);
	end
end
