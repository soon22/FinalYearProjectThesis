clear all;

%37GHz London Plane
a=0.95;
w= 0.95;
Beta = 18;
sigmaT=0.441;
%61.5GHz London Plane
% a=0.25;
% w=0.50;
% Beta = 2 ;
% sigmaT=0.498;

%11 GHz London Plane
% a = 0.70;
% w = 0.95;
% Beta = 100;
% sigmaT = 0.750;

BW3dB=10;
BWr = 0.6*BW3dB;
%Not includes d
Tb=(1-a*w)*sigmaT;
BetaS= 0.6*Beta;
W=(1-a)*w/(1-a*w);

%Characteristic Equation,CE
CE=0;
%N=11
syms s
for n=1:1:10
    CE=CE+(sin(pi./11).*sin(n.*pi./11))./(1-(-cos(n.*pi./11)./s));
end 
n=0;
CE = CE + (sin(pi./22).*sin(pi./22))./(1-(-cos(n.*pi./11)./s));
n=11;
CE = CE + (sin(pi./22)*sin(pi./22))./(1-(-cos(n.*pi./11)./s));
CE = W/2.*CE;
S=solve(CE==1,s);
S=vpa(S);
S=sort(S);
syms A6 A7 A8 A9 A10 A11
Afactor = [A6,A7,A8,A9,A10,A11];

eqn1=0;
for k=6:1:11
    eqn1 = eqn1 + Afactor(k-5) ./ (1-(-cos(6.*pi./11)./S(k+1)));
end
l=eqn1==0;

eqn2=0;
for k=6:1:11
    eqn2 = eqn2 + Afactor(k-5) ./ (1-(-cos(7.*pi./11)./S(k+1)));
end
m=eqn2==0;

eqn3=0;
for k=6:1:11
    eqn3 = eqn3 + Afactor(k-5) ./ (1-(-cos(8.*pi./11)./S(k+1)));
end
n=eqn3==0;

eqn4=0;
for k=6:1:11
    eqn4 = eqn4 + Afactor(k-5) ./ (1-(-cos(9.*pi./11)./S(k+1)));
end
o=eqn4==0;

eqn5=0;
for k=6:1:11
    eqn5 = eqn5 + Afactor(k-5) ./ (1-(-cos(10.*pi./11)./S(k+1)));
end
p=eqn5==0;

eqn6=0;
for k=6:1:11
    eqn6 = eqn6 + Afactor(k-5) ./ (1-(-cos(11.*pi./11)./S(k+1)));
end
q=eqn6==1./(sin(pi./22).*sin(pi./22));

[A,B] = equationsToMatrix([l,m,n,o,p,q], [A6, A7, A8, A9, A10, A11]);
 
Afactors = linsolve(A,B);
Afactors=Afactors.';
sol = solve([l,m,n,o,p,q], [A6, A7, A8, A9, A10, A11]);
Afactor = [sol.A6, sol.A7, sol.A8, sol.A9, sol.A10, sol.A11];
h = 0:0.001:1.5;
r = sqrt(((1.5).^2) - (h.^2));
d = r.*2;
figure(1)
plot(r,h);
figure(2)
plot(d,h);


term1 = exp(-sigmaT.*d);
term2part1 = (exp(-Tb.*d)-exp(-sigmaT.*d)).*4./((BWr.*BWr) + 10.*(BetaS.*BetaS));

term2part2 = 0;
for m = 1:1:10
    term2part2 = term2part2 + (1./factorial(m).*((a.*w.*sigmaT.*d).^m)).*(4./((BWr.*BWr) + m.*(BetaS.*BetaS))-4./((BWr.*BWr) + 10.*(BetaS.*BetaS)));
end
term2part2 = exp(-sigmaT.*d).*term2part2;
term2 = (BWr.*BWr)./4.*(term2part1 + term2part2);
term3part1 = -exp(-Tb.*d).*1./(sin(pi/22).*sin(pi/22));
term3part2 = 0;
for k = 6:1:11
    term3part2 = term3part2 + (Afactors(k-5)).*exp(Tb.*d./S(k+1))./(1-(-cos(pi))./S(k+1));
end
term3 = (BWr.*BWr)./2.*(term3part1 + term3part2);

L = 10.*log10(term1 + term2 + term3);

figure (3)
plot(h,L);