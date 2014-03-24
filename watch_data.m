function watch_data(file_name)
    load(file_name)
    M = 0:0.05:1;
	subplot(2,1,1);
    hist(g_a/g_max,M);
    title('g_a');
	subplot(2,1,2);
	hist(g_a_in/g_max_in,M);
    title('g_a_in');
end