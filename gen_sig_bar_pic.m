function [sig_ex_all, typical] = gen_sig_bar_pic(n,m,tmax,fs,lamda,r_1,tao_c,sigma, IMAGES)
	time_all = fs*tmax;
    
    dim = floor(sqrt(n));

	sig_ex_all = zeros(time_all, n);
    
    lamda_new = lamda*fs;
	tao_c_new = tao_c*fs;
    
    pre_all = random('exp', lamda_new, m+tmax/lamda, round(tmax/lamda*2));
    pre_all = cumsum(pre_all')';
    
    tn = tmax*fs;
    t_left = 1+m;

	for i=1:m
        pre = pre_all(i,:);
        while (pre(end)<tn)
            pre = [pre pre_all(t_left,:)+pre(end)];
            t_left = t_left+1;
        end        
        ans_lists = find((pre>tn),1);
        ans_my = round(pre(1:ans_lists(1)))+1;          
        sig_in_all(ans_my,i) = 1;
	end

    pre_all = random('exp', tao_c_new, 2, round(tmax/tao_c*2));
    pre_sum = cumsum(pre_all')';
    
    t_left = 2;

	pre = pre_sum(1,:);
	interval_my = pre_all(1,:);
	while (pre(end)<tn)
		pre = [pre pre_sum(t_left,:)+pre(end)];
		interval_my = [interval_my pre_all(t_left,:)];
		t_left = t_left+1;
	end        
	ans_lists = find((pre>tn),1);
	interval_my = interval_my(1:ans_lists(1));

    typical     = zeros(dim, dim);
    
	len_inter   = length(interval_my);
    I_0         = 0.4;
    
    sig_dim = dim;
%     sig_dim = 10;
	[x,y,z] = size(IMAGES);
	which_x = ceil(rand(len_inter,1)*(x-sig_dim));
	which_y = ceil(rand(len_inter,1)*(y-sig_dim));
	which_z = ceil(rand(len_inter,1)*z);
    
    r_a = zeros(n, len_inter);
    for j=1:len_inter
        sig_pic_now = IMAGES((which_x(j)):(which_x(j)+sig_dim-1),(which_y(j)):(which_y(j)+sig_dim-1),which_z(j));
        sig_pic_now = sig_pic_now - mean(mean(sig_pic_now));
        sig_pic_now = sig_pic_now/(max(max(max(sig_pic_now)),0.1));
        for i=1:n
            a1  = i-1;
            x1  = floor(a1/dim);
            y1  = mod(a1, dim);
            x2  = floor(x1/(dim/sig_dim));
            y2  = floor(y1/(dim/sig_dim));
            
            fre_now     = sig_pic_now(x2+1, y2+1);
            if fre_now<0.1
                r_a(i,j)    = 0;
            else
%                 r_a(i,j)    = min((fre_now-0.1)/I_0*50, 100);
                r_a(i,j)    = min((fre_now)/I_0*50, 100);
            end

            typical(x1+1, y1+1)     = r_a(i,j);
        end
    end
	r_a = max(r_a, 0.001);

	time_a = repmat(interval_my, n, 1).*r_a/fs;
	max_time = max(max(time_a));
	max_time = ceil(max_time)+1;
	max_time = max_time*2;
	
	lamda_standrad = 0.1*fs;
	big_exp_rand = random('exp', lamda_standrad, n, len_inter, 2*max_time);
	r_a_new = zeros(n,len_inter,2*max_time);
	for i=1:2*max_time
		r_a_new(:,:,i) = r_a;
	end
	big_exp_rand = big_exp_rand.*(fs./r_a_new)/lamda_standrad;

	for i=1:n
		middle_rand = zeros(len_inter, 2*max_time);
		middle_rand(:,:) = big_exp_rand(i,:,:);
		middle_rand = cumsum(middle_rand')';
		new_start = [1 interval_my(1:(end-1))];
		new_start = cumsum(new_start);
		middle_rand = middle_rand + repmat(new_start', 1, 2*max_time);
		new_end = cumsum(interval_my)+1;
		ans_my = find(middle_rand<repmat(new_end',1,2*max_time));
		sig_ex_all(round(middle_rand(ans_my)), i) = 1;
	end
end
