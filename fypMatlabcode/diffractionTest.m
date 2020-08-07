freq = 1e9;
lambda = physconst('LightSpeed')/freq;
d = 100;
h1 = 1/2*t;
h2 = 1/2*t;

a_distances = [h1 , d];
a = norm(a_distances);
b = d;
c_distances = [h2,d ];
c = norm(c_distances);

v1 = h1*sqrt(2/lambda.*(1/a+1/b));
v2 = h2*sqrt(2/lambda.*(1/b+1/c));

C = @(s) cos((pi*(s.^2))/2);
S = @(s) sin((pi*(s.^2))/2);

Cv1 = integral(C,0,v1);
Sv1= integral(S,0,v1);
Cv2 = integral(C,0,v2);
Sv2= integral(S,0,v2);
% Jv is loss
Jv1 = -20.*log(sqrt((1-Cv1-Sv1).^2 + (Cv1-Sv1).^2)/2);
Jv2 = -20.*log(sqrt((1-Cv2-Sv2).^2 + (Cv2-Sv2).^2)/2);
L = Jv1 + Jv2;