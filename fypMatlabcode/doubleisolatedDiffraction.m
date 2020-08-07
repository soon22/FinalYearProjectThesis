%Double isolated diffraction
clear all;
freq = 1e9;
lambda = physconst('LightSpeed')/freq;

% top diffracted
%t - height of tree, d- depth of tree

halft = 50;
a = 100;
b = 100;
c = 100;
topDiffracted = DoubleIsolatedDiffractionClass;
topDiffracted.halft = halft;
topDiffracted.a = a;
topDiffracted.b = b;
topDiffracted.c = c;
topDiffracted.lambda = lambda;
Ltop = calculateLoss(topDiffracted);

w = 100;
%side A diffracted
percentWidth_a = 0.1;
halft_a = percentWidth_a.*w;
side_a_Diffracted = DoubleIsolatedDiffractionClass;
side_a_Diffracted.halft = halft_a;
side_a_Diffracted.a = a;
side_a_Diffracted.b = b;
side_a_Diffracted.c = c;
side_a_Diffracted.lambda = lambda;
Lside_a = calculateLoss(side_a_Diffracted);

%side b diffracted
halft_b = (1-percentWidth_a).*w;
side_b_Diffracted = DoubleIsolatedDiffractionClass;
side_b_Diffracted.halft = halft_b;
side_b_Diffracted.a = a;
side_b_Diffracted.b = b;
side_b_Diffracted.c = c;
side_b_Diffracted.lambda = lambda;
Lside_b = calculateLoss(side_b_Diffracted);





