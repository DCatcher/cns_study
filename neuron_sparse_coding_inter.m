function neuron_sparse_coding_inter(param)
	sig_dim = param.sig_dim;
	ex_n = param.ex_n;
    in_n = param.in_n;
    all_n = ex_n + in_n;
    
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
	sig_n = sig_dim*sig_dim*param.offset;
    vol_ex = zeros(all_n,time_all);
	vol_ex(:,1) = v_reset;
	sta_pic = zeros(all_n, sig_dim, sig_dim);
    sta_pic_ans = zeros(all_n, sig_dim, sig_dim);
	sta_num = zeros(all_n, 1);
    
    if display_mode==1
        sta_figure = figure;
        spike_figure = figure;
        sper_patch_figure = figure;
    end
	for lamda_now = lamda_max
		if (cont==0)
			g_a_sig_to_ex_use = g_max*rand(sig_n, ex_n);
			g_sig_to_ex = zeros(ex_n, 1);
			
% 			g_a_ex_to_other_use = g_max_in*rand(ex_n, ex_n);
% 			g_ex_to_other = 0*ones(ex_n,1);
            g_a_ex_to_in_use = g_max*rand(ex_n, in_n);
            g_ex_to_in = zeros(in_n, 1);
            
            g_a_in_to_other_use = g_max_in*rand(in_n, all_n);
            g_in_to_other = zeros(all_n, 1);
		elseif (cont==1)
			load save_data
        end

%         batch_num_count = 0;
		if time_pattern==1
			pattern_patches = gen_pat_pic_temporal(param);
        end
        
		for patch_now=0:batch_size:train_num
			spike_per_patch = zeros(1, batch_size);

            if short_report_mode==1
                fprintf('patch_now:%i/%i, step: %i\n',patch_now,train_num,batch_size);
            end
            
            if param.one_fire_exp==1
                [sig_ex_all,pic_batch] = gen_sig_one_fire_batch_exp(sig_dim, batch_size, time_per, IMAGES, param.fire_max_time, param.I_0);
            elseif param.one_fire==1
				[sig_ex_all,pic_batch] = gen_sig_one_fire_batch(sig_dim, batch_size, time_per, fs, IMAGES, param.fire_interval);
			elseif time_pattern==0
				[sig_ex_all,pic_batch] = gen_sig_pic(sig_dim,batch_size,time_per,fs,images_in,IMAGES, lamda_max, lamda_min);
			else
				[sig_ex_all,pic_batch] = gen_sig_pic_temporal(sig_dim,batch_size,time_per,fs,images_in, IMAGES, divide_len, pattern_patches);
            end
            
			if (patch_now>0)
				vol_ex(:,1) = vol_ex(:,time_all);
			end
			batch_indx = 1;
			batch_fire_num = zeros(all_n, batch_size);
            g_a_sig_to_ex = g_a_sig_to_ex_use;
            g_a_ex_to_in = g_a_ex_to_in_use;
            g_a_in_to_other = g_a_in_to_other_use;
            
            sig_list = [];
            fire_ex_list = [];
            fire_in_list = [];
			for i=2:time_all
				if mod(i, time_per)==1
					batch_fire_num(:,batch_indx) = sum((vol_ex(:,(i-time_per):(i-1))>-10)',1)';
					spike_per_patch(batch_indx) = sum(batch_fire_num(:,batch_indx)>0);
					for j=1:all_n
						sta_pic(j,:,:) = batch_fire_num(j,batch_indx)*pic_batch(batch_indx,:,:)+sta_pic(j,:,:);
						sta_num(j) = sta_num(j)+batch_fire_num(j,batch_indx);
                    end
					
					vol_ex(:,i-1) = v_reset;
                    batch_indx = batch_indx+1;
                    
                    for fire_his=fire_in_list
                        num_ner_in = fire_his(1);
                        num_time_in = fire_his(2);
                        
                        for fire_his_2 = fire_in_list
                            num_ner_in_2 = fire_his_2(1);
                            num_time_in_2 = fire_his_2(2);
                            
                            if num_ner_in~=num_ner_in_2
                                del_time = num_time_in - num_time_in_2;
                                
                                if del_time<0
                                    g_a_in_to_other(num_ner_in, num_ner_in_2+ex_n) = g_a_in_to_other(num_ner_in, num_ner_in_2+ex_n) + (exp((del_time/fs)/tao_neg_in)-neg_amp)*A_neg_in*g_max_in;
                                else
                                    g_a_in_to_other(num_ner_in, num_ner_in_2+ex_n) = g_a_in_to_other(num_ner_in, num_ner_in_2+ex_n) + (exp((-del_time/fs)/tao_pos_in)-neg_amp)*A_pos_in*g_max_in;
                                end
                            end
                        end
                    end
                    
                    for fire_his=fire_ex_list
                        num_ner_ex = fire_his(1);
                        num_time_ex = fire_his(2);
                        
                        for fire_his_in=fire_in_list
                            num_ner_in = fire_his_in(1);
                            num_time_in = fire_his_in(2);

                            del_time = num_time_ex-num_time_in;
                            
                            if del_time>=0
                                g_a_ex_to_in(num_ner_ex, num_ner_in) = g_a_ex_to_in(num_ner_ex, num_ner_in) - (exp((-del_time/fs)/tao_pos))*A_pos*g_max;
                            else
                                g_a_ex_to_in(num_ner_ex, num_ner_in) = g_a_ex_to_in(num_ner_ex, num_ner_in) + (exp(( del_time/fs)/tao_neg))*A_neg*g_max;
                            end
                            
                            del_time = num_time_in - num_time_ex;

                            if del_time<0
                                g_a_in_to_other(num_ner_in, num_ner_ex) = g_a_in_to_other(num_ner_in, num_ner_ex) + (exp((del_time/fs)/tao_neg_in)-neg_amp)*A_neg_in*g_max_in;
                            else
                                g_a_in_to_other(num_ner_in, num_ner_ex) = g_a_in_to_other(num_ner_in, num_ner_ex) + (exp((-del_time/fs)/tao_pos_in)-neg_amp)*A_pos_in*g_max_in;
                            end
                            
                        end
                        
                        for sig_his_other=sig_list
                            num_ner_sig = sig_his_other(1);
                            num_time_sig = sig_his_other(2);
                            
                            del_time = num_time_sig - num_time_ex;
                            if del_time>=0
                                g_a_sig_to_ex(num_ner_sig, num_ner_ex) = g_a_sig_to_ex(num_ner_sig, num_ner_ex) - (exp((-del_time/fs)/tao_pos))*A_pos*g_max;
                            else
                                g_a_sig_to_ex(num_ner_sig, num_ner_ex) = g_a_sig_to_ex(num_ner_sig, num_ner_ex) + (exp(( del_time/fs)/tao_neg))*A_neg*g_max;
                            end 
                        end
                    end
                    
                    fire_in_list = [];
                    fire_ex_list = [];
                    sig_list = [];
					g_sig_to_ex = zeros(ex_n,1);
					g_ex_to_in = zeros(in_n,1);
                    g_in_to_other = zeros(all_n, 1);
                    
                    g_a_ex_to_in = max(g_a_ex_to_in, 0);
                    g_a_ex_to_in = min(g_a_ex_to_in, g_max);
                    
                    g_a_sig_to_ex = max(g_a_sig_to_ex, 0);
                    g_a_sig_to_ex = min(g_a_sig_to_ex, g_max);
                    
                    g_a_in_to_other = max(g_a_in_to_other, 0);
                    g_a_in_to_other = min(g_a_in_to_other, g_max_in);
                    
                end
                
                vol_ex(:,i) = v_reset*(vol_ex(:,i-1)>v_th) + (vol_ex(:,i-1)<v_th).*...
                    ((v_rest-vol_ex(:,i-1)+[g_sig_to_ex; g_ex_to_in].*(e_ex-vol_ex(:,i-1))+g_in_to_other.*(e_in-vol_ex(:,i-1)))/tao_m*1.0/fs+vol_ex(:,i-1));

				Po_ex = (vol_ex(:,i)>=v_th);
				Po_ex_row = find(Po_ex==1);
				vol_ex(:,i) = vol_ex(:,i).*(vol_ex(:,i)<v_th);
                vol_ex(:,i) = max(e_in, vol_ex(:,i));
                
                
                for Po_ex_tmp = Po_ex_row'
                    if isempty(Po_ex_tmp)
                        break;
                    end
                    if Po_ex_tmp<=ex_n
                        fire_ex_list = [fire_ex_list [Po_ex_tmp;i]];
                    else
                        fire_in_list = [fire_in_list [(Po_ex_tmp-ex_n);i]];
                    end
                end

%sig_to_ex begin
				sig_in = sig_ex_all(i,:);
				sig_in_row = find(sig_in==1);
                for sig_in_tmp = sig_in_row
                    if isempty(sig_in_tmp)
                        break;
                    end
                    sig_list = [sig_list [sig_in_tmp;i]];
                end

				g_sig_to_ex = g_sig_to_ex*(1-1/tao_ex*1.0/fs);
                

				if batch_learn==1
                    g_sig_to_ex = g_sig_to_ex + (sum((g_a_sig_to_ex_use(sig_in_row,:)),1))';
                else
                    g_sig_to_ex = g_sig_to_ex + (sum((g_a_sig_to_ex(sig_in_row,:)),1))';
                end
                
%sig_to_ex finish

%ex_to_other begin
% 				Po_all = Po_ex;

				sig_in = vol_ex(1:ex_n,i)'>v_th;
				sig_in_row = find(sig_in==1);

                g_ex_to_in = g_ex_to_in*(1-1/tao_ex*1.0/fs);
				
                if batch_learn==1
                    g_ex_to_in = g_ex_to_in + (sum((g_a_ex_to_in_use(sig_in_row,:)),1))';
                else
                    g_ex_to_in = g_ex_to_in + (sum((g_a_ex_to_in(sig_in_row,:)),1))';
                end
                
                sig_in = vol_ex((ex_n+1):all_n, i)'>v_th;
                sig_in_row = find(sig_in==1);
                
                g_in_to_other = g_in_to_other*(1-1/(tao_ex_in*fs));
                
                if batch_learn==1
                    g_in_to_other = g_in_to_other + (sum((g_a_in_to_other_use(sig_in_row,:)),1))';
                else
                    g_in_to_other = g_in_to_other + (sum((g_a_in_to_other(sig_in_row,:)),1))';
                end
%in_to_other finish
            end
            
            if batch_learn==1
                g_a_sig_to_ex_use = g_a_sig_to_ex;
                g_a_ex_to_in_use = g_a_ex_to_in;
                g_a_in_to_other_use = g_a_in_to_other;
            end
            
			spike_every = sum(([vol_ex(1:ex_n,:)']>-10));
			ans_my = mean(sum(([vol_ex(1:ex_n,:)']>-10)));
			spike_every_in = sum(([vol_ex((ex_n+1):all_n,:)']>-10));
			ans_my_in = mean(sum(([vol_ex((ex_n+1):all_n,:)']>-10)));            
            if short_report_mode==1
                fprintf('spiking_time:%f, %f\n',ans_my, ans_my_in);
            end
			if (display_mode==1) && (mod(patch_now,1000)==0)
                
				set(0, 'CurrentFigure', sta_figure);
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
				set(0, 'CurrentFigure', spike_figure);
				hist(spike_every);
				set(0, 'CurrentFigure', sper_patch_figure);
				hist(spike_per_patch);
                
				pause(1);
            end
            fprintf('Zero Patch No: %f\n',sum(spike_per_patch==0));
            
            tmp_sta = find(sta_num>param.sta_big_num);
            for tmp_idx=tmp_sta'
                sta_pic_ans(tmp_idx,:,:) = sta_pic(tmp_idx,:,:)./sta_num(tmp_idx);
            end
		end
		sig_ex_all = [];
        s = ['data_',int2str(ex_n),'_',int2str(in_n)];
        save(s);
        fprintf('saved!%s\n',s);
	end
end
