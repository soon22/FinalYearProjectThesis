close all;
clear all;
psi = 0.01:0.05:90;
[rh,rv] = ref_coef(psi,65,30.7);
gamamodv = abs(rv);
gamamodh = abs(rh);
subplot(2,1,1)
plot(psi,gamamodv,'k',psi,gamamodh,'k-','linewidth',1.5);
grid
legend('Vertical Polarization','Horizontal Polarization')
title('Reflection coefficient - magnitude')
pv = -angle(rv);
ph = angle(rh);
subplot(2,1,2)
plot(psi,pv,'k',psi,ph,'k-','linewidth',1.5);
grid
legend('Vertical Polarization','Vertical Polarization')
title('Reflection coefficient -phase')
xlabel('Grazing angle in degrees')

