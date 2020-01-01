

fsamp = 2800*2*pi;
fcuts = [500*2*pi 600*2*pi 900*2*pi 1050*2*pi];
mags = [1 0 1];
devs = [0.05 0.0044660 0.05];

[n,Wn,beta,ftype] = kaiserord(fcuts,mags,devs,fsamp);
hh = fir1(n,Wn,ftype,kaiser(n+1,beta),'noscale');

freqz(hh,1,1024,fsamp)