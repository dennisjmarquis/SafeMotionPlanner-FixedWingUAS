% load ekfTestF9.mat
load F3_nofilt.mat

Fs = 1/dt;
Fc = 10;
[den, num] = butter(4, Fc/(Fs/2), 'low');  % 4th-order low-pass Butterworth filter


R_p =var(true.p-filtfilt(den,num,true.p))
R_q =var(true.q-filtfilt(den,num,true.q))
R_r =var(true.r-filtfilt(den,num,true.r))

R_u =var(true.u-filtfilt(den,num,true.u))
R_v =var(true.v-filtfilt(den,num,true.v))
R_w =var(true.w-filtfilt(den,num,true.w))

R_phi =var(true.phi-filtfilt(den,num,true.phi))
R_theta =var(true.theta-filtfilt(den,num,true.theta))
R_psi =var(true.psi-filtfilt(den,num,true.psi))

R_va =var(true.Va-filtfilt(den,num,true.Va))

R_ax =var(true.ax-filtfilt(den,num,true.ax))
R_ay =var(true.ay-filtfilt(den,num,true.ay))
R_az =var(true.az-filtfilt(den,num,true.az))
