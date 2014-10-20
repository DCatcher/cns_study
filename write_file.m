function param   = write_file(param)

if mod(param.time_now, param.save_per_time_m)~=0
    return;
end

param.wirten        = 1;

param.speed_ex      = [];
param.vol_ex        = [];

save(param.save_file_name, 'param');