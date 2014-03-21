initial

n = 100;
m = 20;
tmax = 5;
fs = 1000;

[vol, g_a] = neuron(n,m,tmax,fs,0);

M = 0:0.05:1;
%hist(g_a/0.015,M);
%title('g_a');
ans = sum(vol>-10);
disp(ans);
