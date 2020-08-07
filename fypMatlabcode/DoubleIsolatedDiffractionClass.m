classdef DoubleIsolatedDiffractionClass
    properties
        a;
        b;
        c;
        halft;
        lambda;

    end
    methods
        function L = calculateLoss(obj)
            h1 = obj.halft - obj.halft*obj.a/(obj.a+obj.b);
            h2 = obj.halft - obj.halft*obj.c/(obj.c+obj.b);
            d2 = obj.b;
            d1 = norm([obj.a,obj.halft]);
            d3 = norm([obj.c,obj.halft]);
            v1 = h1*sqrt(2/obj.lambda*(1/d1+1/d2));
            v2 = h2*sqrt(2/obj.lambda*(1/d2+1/d3));
             
            C = @(s) cos((pi*(s.^2))/2);
            S = @(s) sin((pi*(s.^2))/2);
            
            Cv1 = integral(C,0,v1);
            Sv1= integral(S,0,v1);
            Cv2 = integral(C,0,v2);
            Sv2= integral(S,0,v2);
            % Jv is loss
            Jv1 = -20.*log(sqrt((1-Cv1-Sv1).^2 + (Cv1-Sv1).^2)/2);
            Jv2 = -20.*log(sqrt((1-Cv2-Sv2).^2 + (Cv2-Sv2).^2)/2);
            Lc = 10.*log((obj.a+obj.b)*(obj.b+obj.c)/(obj.b*(obj.a+obj.b+obj.c)));
            
            L = Jv1 + Jv2 + Lc;
        end
    end
end

