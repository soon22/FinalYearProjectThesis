a=0.95;
w=0.95;
Beta = 18;
sigmaT=0.441;
BW3dB=1.2;
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
    CE=CE+(sin(pi/11)*sin(n*pi/11))/(1-(cos(n*pi/11)/s));
end 
n=0;
CE = CE + (sin(pi/22)*sin(pi/22))/(1-(cos(n*pi/11)/s));
n=11;
CE = CE + (sin(pi/22)*sin(pi/22))/(1-(cos(n*pi/11)/s));
CE = W/2*CE;
S=solve(CE==1,s);
S=vpa(S);
S=sort(S);
syms A6 A7 A8 A9 A10 A11
Afactor = [A6,A7,A8,A9,A10,A11];

eqn1=0;
for k=6:1:11
    eqn1 = eqn1 + Afactor(k-5) / (1-(-cos(6*pi/11)/S(k+1)));
end
l=eqn1==0;

eqn2=0;
for k=6:1:11
    eqn2 = eqn2 + Afactor(k-5) / (1-(-cos(7*pi/11)/S(k+1)));
end
m=eqn2==0;

eqn3=0;
for k=6:1:11
    eqn3 = eqn3 + Afactor(k-5) / (1-(-cos(8*pi/11)/S(k+1)));
end
n=eqn3==0;

eqn4=0;
for k=6:1:11
    eqn4 = eqn4 + Afactor(k-5) / (1-(-cos(9*pi/11)/S(k+1)));
end
o=eqn4==0;

eqn5=0;
for k=6:1:11
    eqn5 = eqn5 + Afactor(k-5) / (1-(-cos(10*pi/11)/S(k+1)));
end
p=eqn5==0;

eqn6=0;
for k=6:1:11
    eqn6 = eqn6 + Afactor(k-5) / (1-(-cos(11*pi/11)/S(k+1)));
end
q=eqn6==1/(sin(pi/22)*sin(pi/22));

[A,B] = equationsToMatrix([l,m,n,o,p,q], [A6, A7, A8, A9, A10, A11]);
 
Afactor = linsolve(A,B);
Afactor=Afactor.';
d=0:100;
Iri = exp(-sigmaT*d);
I1term1 = (exp(-Tb*d)-exp(-sigmaT*d))*4/(BWr*BWr + 10*BetaS*BetaS);
I1term2=0;
for m=1:1:10
    I1term2=I1term2+1/factorial(m)*((a*w*sigmaT*d).^m)*(4/(BWr*BWr+m*BetaS*BetaS)-4/(BWr*BWr+10*BetaS*BetaS));
end
I1=BWr*BWr/4*(I1term1+(exp(-sigmaT*d)).*I1term2);

I2term1=-exp(Tb*d)*1/(sin(pi/22)*sin(pi/22));
I2term2=0;
for k=6:1:11
   I2term2=I2term2+ (Afactor(k-5))*exp(Tb*d/S(k+1))*1/(1-(-cos(11/11*pi))/S(k+1));
end
I2=BWr*BWr/2*(I2term1+I2term2);
L=-10*log10(Iri+I1+I2);


plot(d,L);
