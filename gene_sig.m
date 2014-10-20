function param   = gene_sig(param)

if param.with_speed==0
    disp('Not with speed not finished');
    assert(false);
end

if param.with_speed==1
    tmp_file_name   = [param.save_sig_file, '/', param.save_sig_name, '_', int2str(param.time_now)];
    
    if param.sig_from_file==0
        if param.choose_rand==0
            [tmp_sig_ex, tmp_sig_in, speed_ex] = gen_sig_ex_with_speed(param.sig_n,0,param.tmax,param.fs,param.lamda);
        else
            test_i                      = randsample(1:size(param.lamda, 2), 1);
            param.rand_choices(end+1)   = test_i;
            
            [tmp_sig_ex, tmp_sig_in, speed_ex] = gen_sig_ex_with_speed(param.sig_n,0,param.tmax,param.fs,param.lamda(:, test_i));
        end
    else
        load(tmp_file_name);
    end
    
    if param.save_sig==1
        save(tmp_file_name, 'speed_ex');
    end
    
    param.speed_ex      = speed_ex;
end