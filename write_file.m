function param   = write_file(param)

if mod(param.time_now, param.save_per_time)~=0
    return;
end

param.speed_ex      = [];
param.vol_ex        = [];

save(param.save_file_name, 'param');