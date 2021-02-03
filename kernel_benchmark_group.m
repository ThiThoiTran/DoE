function [f] = kernel_benchmark_group(x,n_level,i_level,mode)
% function to create the benchmark function (~ kernel herding paper) but
% for necklace prob
x = x(:);
NC = size(x,1);
if mode == 1
    xl = 0;
    xu = 5;
    X = zeros(n_level,1);
    for i = 1:n_level
        X(i) = xl + (i-1)*(xu - xl)/(n_level-1);
    end
    f = X(i_level)*ones(NC,1) + x;
elseif mode == 2
    xl = 0;
    xu = 5;
    X = zeros(n_level,1);
    for i = 1:n_level
        X(i) = xl + (i-1)*(xu - xl)/(n_level-1);
    end
    f = (X(i_level)*ones(NC,1) + x).^2;
elseif mode == 3
    xl = 0;
    xu = 5;
    X = zeros(n_level,1);
    for i = 1:n_level
        X(i) = xl + (i-1)*(xu - xl)/(n_level-1);
    end
    f = (X(i_level)*ones(NC,1) + x).^3;
elseif mode == 4
    xl = 0;
    xu = 5;
    X = zeros(n_level,1);
    for i = 1:n_level
        X(i) = xl + (i-1)*(xu - xl)/(n_level-1);
    end
    f =  sin(norm(x+X(i_level)*ones(NC,1)));
elseif mode == 5 % Branin
    xl = -5;
    xu = 10;
    X = zeros(n_level,1);
    for i = 1:n_level
        X(i) = xl + (i-1)*(xu - xl)/(n_level-1);
    end
    f = (X(i_level)-5.1/(4*pi^2)*x^2+5/(pi)*x-6)^2+10*(1-1/(8*pi))*cos(x)+10;
elseif mode == 6 % hartman3
    xl = 0;
    xu = 1;
    X = zeros(n_level,1);
    for i = 1:n_level
        X(i) = xl + (i-1)*(xu - xl)/(n_level-1);
    end
    a = [ [3.0,  0.1,  3.0,  0.1];
        [10.0, 10.0, 10.0, 10.0];
        [30.0, 35.0, 30.0, 35.0] ];
    p = [ [0.36890, 0.46990, 0.10910, 0.03815];
        [0.11700, 0.43870, 0.87320, 0.57430];
        [0.26730, 0.74700, 0.55470, 0.88280] ];
    c = [1.0, 1.2, 3.0, 3.2];
    x(3) = X(i_level);
    T = 0;
    for k = 1:4
        K = 0;
        for j = 1:3
            K = K + a(j,k)*(x(j)-p(j,k))^2;
        end
        T = T +c(k)*exp(-K);
    end
    f = -T;
elseif mode == 7 %perm8
    dimen = 8;
    xl = 0;
    xu = 3;
    beta = 100;
    X = zeros(n_level,1);
    for i = 1:n_level
        X(i) = xl + (i-1)*(xu - xl)/(n_level-1);
    end
    x(dimen) = X(i_level);
    T = 0;
    for k = 1:8
        S = 0;
        for j = 1:8
            S = S + ((j+1)+beta)*(x(j)^k -(1/(j+1))^k);
        end
        T = T + S^2+1000;
    end
    f = T; 
elseif mode == 8 % Hartman6
     dimen = 6;
    xl = 0;
    xu = 3;
    X = zeros(n_level,1);
    for i = 1:n_level
        X(i) = xl + (i-1)*(xu - xl)/(n_level-1);
    end
    a = [ [10.00,  0.05,  3.00, 17.00];
        [3.00, 10.00,  3.50,  8.00];
        [17.00, 17.00,  1.70,  0.05];
        [3.50,  0.10, 10.00, 10.00];
        [1.70,  8.00, 17.00,  0.10];
        [8.00, 14.00,  8.00, 14.00] ];
    p = [ [0.1312, 0.2329, 0.2348, 0.4047];
        [0.1696, 0.4135, 0.1451, 0.8828];
        [0.5569, 0.8307, 0.3522, 0.8732];
        [0.0124, 0.3736, 0.2883, 0.5743];
        [0.8283, 0.1004, 0.3047, 0.1091];
        [0.5886, 0.9991, 0.6650, 0.0381] ];
    c = [1.0, 1.2, 3.0, 3.2];
    
        x(dimen) = X(i_level);
        T = 0;
        for i = 1:4
            K = 0;
            for j = 1:6
                K = K + a(j,i)*(x(j)-p(j,i))^2;
            end
            T = T +c(i)*exp(-K);
        end
        f = -T;
elseif mode == 9 % perm6
        dimen = 6;
        xl = -6;
        xu = 6;
        beta = 60;
        X = zeros(n_level,1);
        for i = 1:n_level
            X(i) = xl + (i-1)*(xu - xl)/(n_level-1);
        end
        x(dimen) = X(i_level);
        T = 0;
        for k = 1:6
            S = 0;
            for i = 1:6
                S = S + ((i+1)^k+beta)*(x(i)/(i+1)^k-1);
            end
            T = T + S^2+1000;
        end
        f = T;
elseif mode == 10 %Wong 2
     ref = x(1)^2+x(2)^2+x(1)*x(2)-14*x(1)-16*x(2)+(x(3)-10)^2+4*(x(4)-5)^2+(x(5)-3)^2 ...
        +2*(x(6)-1)^2 + 5*x(7)^2+7*(x(8)-11)^2+2*(x(9)-10)^2+(x(10)-7)^2+45;
       if i_level == 1 
           f = ref;
       elseif i_level == 2 
           f = ref + 10*(3*(x(1)-2)^2+4*(x(2)-3)^2+2*x(3)^2-7*x(4)-120);
       elseif i_level == 3
           f = ref + 10*(5*(x(1))^2+8*(x(2))^2+(x(3)-6)^2-2*x(4)-40);
       elseif i_level == 4 
          f = ref + 10*(0.5*(x(1)-8)^2+2*(x(2)-4)^2+3*x(5)^2-x(6)-30);
       elseif i_level == 5 
            f = ref + 10*((x(1))^2+2*(x(2)-2)^2-2*x(1)*x(2)+14*x(5)-6*x(6));
       elseif i_level == 6
           f = ref + 10*(-3*(x(1))+6*(x(2))+12*(x(9)-8)^2-7*x(10));
       end
elseif mode == 11 %Wong3
   %  M = countSol(6);
       ref = x(1)^2+x(2)^2+x(1)*x(2)-14*x(1)-16*x(2)+(x(3)-10)^2+4*(x(4)-5)^2+(x(5)-3)^2 ...
        +2*(x(6)-1)^2 + 5*x(7)^2+7*(x(8)-11)^2+2*(x(9)-10)^2+(x(10)-7)^2+(x(11)-9)^2 ...
        +10*(x(12)-1)^2+5*(x(13)-7)^2+4*(x(14)-14)^2+27*(x(15)-1)^2+x(16)^4 ...
        +(x(17)-2)^2+13*(x(18)-2)^2+(x(19)-3)^2+x(20)^2+95;
        if i_level == 1 
            f = ref;
        elseif i_level == 2 
            f = ref+10*(3*(x(1)-2)^2+4*(x(2)-3)^2+2*x(3)^2-7*x(4)-120);
        elseif i_level == 3 
            f = ref + 10*(5*(x(1))^2+8*(x(2))^2+(x(3)-6)^2-2*x(4)-40);
        elseif i_level == 4 
            f = ref + 10*(0.5*(x(1)-8)^2+2*(x(2)-4)^2+3*x(5)^2-x(6)-30);
       elseif i_level == 5 
            f = ref + 10*((x(1))^2+2*(x(2)-2)^2-2*x(1)*x(2)+14*x(5)-6*x(6));
        elseif i_level == 6 
            f = ref + 10*(-3*(x(1))+6*(x(2))+12*(x(9)-8)^2-7*x(10));
       elseif i_level == 7 
            f = ref + 10*(x(1)^2+5*x(11)-8*x(12)-28);
       elseif i_level == 8 
            f = ref + 10*(4*x(1)+9*x(2)+5*x(13)^2-9*x(14)-87);
       elseif i_level == 9 
            f = ref + 10*(3*(x(1))+4*(x(2))+3*(x(13)-6)^2-14*x(14)-10);
       elseif i_level == 10 
            f = ref + 10*(10*x(1)^2+35*x(15)-79*x(16)-92);
       elseif i_level == 11 
            f = ref + 10*(15*x(2)^2+11*x(15)-61*x(16)-54);
       elseif i_level == 12 
            f = ref + 10*(5*x(1)^2+2*x(2)+9*x(17)^4-x(18)-68);
       elseif i_level == 13 
            f = ref + 10*(x(1)^2-x(9)+19*x(19)-20*x(20)+19);
       elseif i_level == 14 
            f = ref + 10*(7*x(1)^2+5*x(2)^2+x(19)^2-30*x(20));
        end
end
end


