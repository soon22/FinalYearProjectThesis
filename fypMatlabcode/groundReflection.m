clear all;
%Calculate reflection coefficient
%relative permittivity and conductivity 

%dry ground
% relPerm = 5;
% cond = 10;

%medium dry ground
%relPerm = 5;
%cond = 8;

%  wet ground
relPerm = 3;
cond = 0.45;


freq = 37e9;

lambda = physconst('LightSpeed')/freq;
psi = 0.01:0.05:90;
angle = psi.*pi./180;
%calculate complex permittivity
complexPerm = relPerm - 1i*60*lambda*cond;

% for horizontal polarization
Ch = complexPerm - (cos(angle)).^2;
%for vertical polarization
Cv = (complexPerm - (cos(angle)).^2)./((complexPerm).^2);

RoVer = (sin(angle) - sqrt(Cv))./(sin(angle) + sqrt(Cv));
RoHor = (sin(angle) - sqrt(Ch))./(sin(angle) + sqrt(Ch));

brewsterAngle = radtodeg(asin(1/sqrt(abs(complexPerm))));
%plot(psi,abs(RoVer),psi,abs(RoHor));
d0 = 103;
d1 = sin(angle).*103./sin(pi-2.*angle);
d2 = d1;
Lgroundver = 20*log((d1+d2)/d0)-20*log(abs(RoVer));

Lgroundhor = 20*log((d1+d2)/d0)-20*log(abs(RoHor));
figure(1)
plot(psi,Lgroundver)
figure(2)
plot(psi,Lgroundhor)




