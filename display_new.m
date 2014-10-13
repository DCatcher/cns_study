function param = display_new(param, screen_new)

if nargin < 2
    screen_new  = 0;
end

if mod(param.time_now, param.display_per_time)~=0
    return;
end

if param.display_mode==0
    return;
end

if (param.screen_initial==0) || (screen_new==1)
    param.screen_initial    = 1;
    param.screen_figure     = figure;
end

set(0, 'CurrentFigure', param.screen_figure);

plot_x  = floor(sqrt(param.neuron_n));
while mod(param.neuron_n, plot_x)~=0
    plot_x  = plot_x - 1;
end
plot_y  = param.neuron_n/plot_x;

for i=1:param.neuron_n
    g_tmp       = param.g_a_sig_to_ex(i, :)/param.g_max;
    g_tmp_new   = zeros(1, length(g_tmp)/param.display_cond);
    for j=1:length(g_tmp)/param.display_cond
        g_tmp_new(j)    = sum(g_tmp((j-1)*param.display_cond+1:j*param.display_cond));
    end
    subplot(plot_x, plot_y, i);
    plot(g_tmp_new);
end

pause(0.01);