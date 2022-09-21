close all
clear
clc

% path and name of txt int count files
V=load('1diff_0_1.txt');
F=load('1diff_45_1.txt');
H=load('1diff_90_1.txt');
T=load('1diff_135_1.txt');

stWL=380; % starting wavelength
endWL=700; % end wavelength

W=V(:,1);
w=W(W>380&W<700); % selecged wavelength region

v=V(W>380&W<700,2); % vertical (0)
f=F(W>380&W<700,2); % 45
h=H(W>380&W<700,2); % horizontal (90)
t=T(W>380&W<700,2); % 135

S1=(v-h)./(v+h); 
S2=(f-t)./(f+t);

D=sqrt((S1.*S1)+(S2.*S2)); % degree across wavelength range
d=nanmean(D); % average degree

A = (180/pi)*(0.5.*atan(S2./S1)); %
a=nanmean(A);

figure, plot(w,A)
title('Angle')
xlabel('wavlength/nm')
ylabel('angle of pol / degrees')

figure, plot(w,D)
title('Degree of Linear Polarization')
xlabel('wavlength/nm')
ylabel('degree')
ylim([0 1])

{'degree of polarization',d}
{'angle of polarization',a}
