function [sig_ex_all, typical] = gen_sig_bar_2d(n,m,tmax,fs,lamda,r_1,tao_c,sigma, gamma)
	time_all = fs*tmax;
    
    dim = floor(sqrt(n));

	sig_ex_all = zeros(time_all, n);
    
    test_speed = 0;
    
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
    
	len_inter = length(interval_my);
	s = rand(1, len_inter)*n;
    
%     t = rand(1, len_inter)*2*pi;
    t = rand(1, len_inter)*pi;
    
    s_2 = rand(1, len_inter)*n;
    
%     t = rand(1, len_inter)*2*pi;
    t_2 = rand(1, len_inter)*pi;

	r_ba = 1/lamda;
	a = (1:n)';
    r_a = zeros(n, len_inter);
% 	r_a = r_ba + r_1*(exp(-(s-a).^2/(2*sigma^2))+exp(-(s-a+n).^2/(2*sigma^2))+exp(-(s-a-n).^2/(2*sigma^2)));
    for i=1:n
        for j=1:len_inter
            a1  = i;
            a2  = s(j);
            x1  = floor(a1/dim);
            y1  = mod(a1, dim);
            x2  = floor(a2/dim);
            y2  = mod(a2, dim);
            
            theta   = t(j);
%             theta   = pi/4;
            dist    = 0;
            
            sigma_x     = sigma;
            sigma_y     = sigma/gamma;
            for l=-1:1
                for k=-1:1
                    x_tmp       =(x1 - x2 - l*dim);
                    y_tmp       =(y1 - y2 - k*dim);
                    x_theta     = x_tmp*cos(theta)+y_tmp*sin(theta);
                    y_theta     =-x_tmp*sin(theta)+y_tmp*cos(theta);
                    
                    dist        = dist + min(5*exp(-.5*(x_theta.^2/sigma_x^2+y_theta.^2/sigma_y^2)), 1);
                end
            end
            
            r_a(i,j)    = r_ba + r_1*dist;
%             
%             a2  = s_2(j);
%             x2  = floor(a2/dim);
%             y2  = mod(a2, dim);
%             
%             theta   = t_2(j);
% %             theta   = pi/4;
%             dist    = 0;
%             
%             sigma_x     = sigma;
%             sigma_y     = sigma/gamma;
%             for l=-1:1
%                 for k=-1:1
%                     x_tmp       =(x1 - x2 - l*dim);
%                     y_tmp       =(y1 - y2 - k*dim);
%                     x_theta     = x_tmp*cos(theta)+y_tmp*sin(theta);
%                     y_theta     =-x_tmp*sin(theta)+y_tmp*cos(theta);
%                     
%                     dist        = dist + min(5*exp(-.5*(x_theta.^2/sigma_x^2+y_theta.^2/sigma_y^2)), 1);
%                 end
%             end
%             
%             r_a(i,j)    = max(r_ba + r_1*dist, r_a(i,j));
            
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
