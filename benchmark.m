function [f] = benchmark(mode,x,y,Neck) % benhmark functions checking DOE (15 benchmarks with 4 extended up to 12 binaries vars)
if mode == 1 % CB2 (2,2)
    if sum(y) == 0
        f = x(1)^2 + x(2)^4;
    elseif sum(y) == 1
        f = (2 - x(1))^2 + (2 - x(2))^2;
    elseif sum(y) == 2
        f = 2*exp(x(2) - x(1));
    end
elseif mode == 2 % Rosen
    if sum(y) == 0
        f = x(1)^2 + x(2)^2 + 2*x(3)^2 + x(4)^2 - 5*x(1) - 5*x(2) - 21*x(3) + 7*x(4);
    elseif sum(y) == 1
        f = x(1)^2 + x(2)^2 + 2*x(3)^2 + x(4)^2 - 5*x(1) - 5*x(2) - 21*x(3) + 7*x(4)...
            + 10*(x(1)^2 + x(2)^2 + x(3)^2 + x(4)^2 + x(1) -x(2) + x(3) -x(4) -8);
    elseif sum(y) == 2
        f = x(1)^2 + x(2)^2 + 2*x(3)^2 + x(4)^2 - 5*x(1) - 5*x(2) - 21*x(3) + 7*x(4)...
            + 10*(x(1)^2 + 2*x(2)^2 + x(3)^2 + 2*x(4)^2 - x(1) -x(4)- 10);
    elseif sum(y) == 3
        f = x(1)^2 + x(2)^2 + 2*x(3)^2 + x(4)^2 - 5*x(1) - 5*x(2) - 21*x(3) + 7*x(4)...
            + 10*(2*x(1)^2 + x(2)^2 + x(3)^2 + 2*x(1) -x(2) -x(4)- 5);
    end
elseif mode == 3 % Pentagon
    if sum(y) == 0
        f = -sqrt((x(1) - x(3))^2 + (x(2) - x(4))^2);
    elseif sum(y) == 1
        f = -sqrt((x(3) - x(5))^2 + (x(4) - x(6))^2);
    elseif sum(y) == 2
        f = -sqrt((x(5) - x(1))^2 + (x(6) - x(2))^2);
    end
elseif mode == 4 %Wong2
    ref = x(1)^2+x(2)^2+x(1)*x(2)-14*x(1)-16*x(2)+(x(3)-10)^2+4*(x(4)-5)^2+(x(5)-3)^2 ...
        +2*(x(6)-1)^2 + 5*x(7)^2+7*(x(8)-11)^2+2*(x(9)-10)^2+(x(10)-7)^2+45;
    if isequal(y,[0 0 0 0])
        f = ref;
    elseif isequal(y,[1 0 0 0]) || isequal(y,[0 1 0 0])|| isequal(y,[0 0 1 0])...
            || isequal(y,[0 0 0 1])
        f = ref + 10*(3*(x(1)-2)^2+4*(x(2)-3)^2+2*x(3)^2-7*x(4)-120);
    elseif isequal(y,[1 0 1 0])|| isequal(y,[0 1 0 1])
        f = ref + 10*(5*(x(1))^2+8*(x(2))^2+(x(3)-6)^2-2*x(4)-40);
    elseif isequal(y,[1 1 0 0])||isequal(y,[0 1 1 0])||isequal(y,[0 0 1 1])...
            ||isequal(y,[1 0 0 1])
        f = ref + 10*(0.5*(x(1)-8)^2+2*(x(2)-4)^2+3*x(5)^2-x(6)-30);
    elseif isequal(y,[1 1 1 0])||isequal(y,[0 1 1 1])||isequal(y,[1 0 1 1])...
            ||isequal(y,[1 1 0 1])
        f = ref + 10*((x(1))^2+2*(x(2)-2)^2-2*x(1)*x(2)+14*x(5)-6*x(6));
    elseif isequal(y,[1 1 1 1])
        f = ref + 10*(-3*(x(1))+6*(x(2))+12*(x(9)-8)^2-7*x(10));
    end
elseif mode == 5 % Wong3 (20,6)
       M = countSol(6);
       ref = x(1)^2+x(2)^2+x(1)*x(2)-14*x(1)-16*x(2)+(x(3)-10)^2+4*(x(4)-5)^2+(x(5)-3)^2 ...
        +2*(x(6)-1)^2 + 5*x(7)^2+7*(x(8)-11)^2+2*(x(9)-10)^2+(x(10)-7)^2+(x(11)-9)^2 ...
        +10*(x(12)-1)^2+5*(x(13)-7)^2+4*(x(14)-14)^2+27*(x(15)-1)^2+x(16)^4 ...
        +(x(17)-2)^2+13*(x(18)-2)^2+(x(19)-3)^2+x(20)^2+95;
        if isequal(y,M(1,:))
            f = ref;
        elseif isequal(y,M(2,:))||isequal(y,circshift(M(2,:),[0,1]))||isequal(y,circshift(M(2,:),[0,2]))||isequal(y,circshift(M(2,:),[0,3])) ...
                || isequal(y,circshift(M(2,:),[0,4]))|| isequal(y,circshift(M(2,:),[0,5]))||isequal(y,circshift(M(2,:),[0,6]))
            f = ref+10*(3*(x(1)-2)^2+4*(x(2)-3)^2+2*x(3)^2-7*x(4)-120);
        elseif isequal(y,M(3,:))||isequal(y,circshift(M(3,:),[0,1]))||isequal(y,circshift(M(3,:),[0,2]))||isequal(y,circshift(M(3,:),[0,3])) ...
                || isequal(y,circshift(M(3,:),[0,4]))|| isequal(y,circshift(M(3,:),[0,5]))||isequal(y,circshift(M(3,:),[0,6]))
            f = ref + 10*(5*(x(1))^2+8*(x(2))^2+(x(3)-6)^2-2*x(4)-40);
        elseif isequal(y,M(4,:))||isequal(y,circshift(M(4,:),[0,1]))||isequal(y,circshift(M(4,:),[0,2]))||isequal(y,circshift(M(4,:),[0,3])) ...
                || isequal(y,circshift(M(4,:),[0,4]))|| isequal(y,circshift(M(4,:),[0,5]))||isequal(y,circshift(M(4,:),[0,6]))
            f = ref + 10*(0.5*(x(1)-8)^2+2*(x(2)-4)^2+3*x(5)^2-x(6)-30);
        elseif isequal(y,M(5,:))||isequal(y,circshift(M(5,:),[0,1]))||isequal(y,circshift(M(5,:),[0,2]))||isequal(y,circshift(M(5,:),[0,3])) ...
                || isequal(y,circshift(M(5,:),[0,4]))|| isequal(y,circshift(M(5,:),[0,5]))||isequal(y,circshift(M(5,:),[0,6]))
            f = ref + 10*((x(1))^2+2*(x(2)-2)^2-2*x(1)*x(2)+14*x(5)-6*x(6));
        elseif isequal(y,M(6,:))||isequal(y,circshift(M(6,:),[0,1]))||isequal(y,circshift(M(6,:),[0,2]))||isequal(y,circshift(M(6,:),[0,3])) ...
                || isequal(y,circshift(M(6,:),[0,4]))|| isequal(y,circshift(M(6,:),[0,5]))||isequal(y,circshift(M(6,:),[0,6]))
            f = ref + 10*(-3*(x(1))+6*(x(2))+12*(x(9)-8)^2-7*x(10));
        elseif isequal(y,M(7,:))||isequal(y,circshift(M(7,:),[0,1]))||isequal(y,circshift(M(7,:),[0,2]))||isequal(y,circshift(M(7,:),[0,3])) ...
                || isequal(y,circshift(M(7,:),[0,4]))|| isequal(y,circshift(M(7,:),[0,5]))||isequal(y,circshift(M(7,:),[0,6]))
            f = ref + 10*(x(1)^2+5*x(11)-8*x(12)-28);
        elseif isequal(y,M(8,:))||isequal(y,circshift(M(8,:),[0,1]))||isequal(y,circshift(M(8,:),[0,2]))||isequal(y,circshift(M(8,:),[0,3])) ...
                || isequal(y,circshift(M(8,:),[0,4]))|| isequal(y,circshift(M(8,:),[0,5]))||isequal(y,circshift(M(8,:),[0,6]))

            f = ref + 10*(4*x(1)+9*x(2)+5*x(13)^2-9*x(14)-87);
        elseif isequal(y,M(9,:))||isequal(y,circshift(M(9,:),[0,1]))||isequal(y,circshift(M(9,:),[0,2]))||isequal(y,circshift(M(9,:),[0,3])) ...
                || isequal(y,circshift(M(9,:),[0,4]))|| isequal(y,circshift(M(9,:),[0,5]))||isequal(y,circshift(M(9,:),[0,6]))
            f = ref + 10*(3*(x(1))+4*(x(2))+3*(x(13)-6)^2-14*x(14)-10);
        elseif isequal(y,M(10,:))||isequal(y,circshift(M(10,:),[0,1]))||isequal(y,circshift(M(10,:),[0,2]))||isequal(y,circshift(M(10,:),[0,3])) ...
                || isequal(y,circshift(M(10,:),[0,4]))|| isequal(y,circshift(M(10,:),[0,5]))||isequal(y,circshift(M(10,:),[0,6]))
            f = ref + 10*(10*x(1)^2+35*x(15)-79*x(16)-92);
        elseif isequal(y,M(11,:))||isequal(y,circshift(M(11,:),[0,1]))||isequal(y,circshift(M(11,:),[0,2]))||isequal(y,circshift(M(11,:),[0,3])) ...
                || isequal(y,circshift(M(11,:),[0,4]))|| isequal(y,circshift(M(11,:),[0,5]))||isequal(y,circshift(M(11,:),[0,6]))
            f = ref + 10*(15*x(2)^2+11*x(15)-61*x(16)-54);
        elseif isequal(y,M(12,:))||isequal(y,circshift(M(12,:),[0,1]))||isequal(y,circshift(M(12,:),[0,2]))||isequal(y,circshift(M(12,:),[0,3])) ...
                || isequal(y,circshift(M(12,:),[0,4]))|| isequal(y,circshift(M(12,:),[0,5]))||isequal(y,circshift(M(12,:),[0,6]))

            f = ref + 10*(5*x(1)^2+2*x(2)+9*x(17)^4-x(18)-68);
        elseif isequal(y,M(13,:))||isequal(y,circshift(M(13,:),[0,1]))||isequal(y,circshift(M(13,:),[0,2]))||isequal(y,circshift(M(13,:),[0,3])) ...
                || isequal(y,circshift(M(13,:),[0,4]))|| isequal(y,circshift(M(13,:),[0,5]))||isequal(y,circshift(M(13,:),[0,6]))
            f = ref + 10*(x(1)^2-x(9)+19*x(19)-20*x(20)+19);
        elseif isequal(y,M(14,:))
            f = ref + 10*(7*x(1)^2+5*x(2)^2+x(19)^2-30*x(20));
        end
elseif mode == 6 % Branin
    n_level = 4;
    xl = -5;
    xu = 10;
    X = zeros(n_level,1);
    for i = 1:n_level
        X(i) = xl + (i-1)*(xu - xl)/(n_level-1);
    end
    if sum(y) == 0
        f = (X(1)-5.1/(4*pi^2)*x^2+5/(pi)*x-6)^2+10*(1-1/(8*pi))*cos(x)+10;
    elseif sum(y) == 1
        f = (X(2)-5.1/(4*pi^2)*x^2+5/(pi)*x-6)^2+10*(1-1/(8*pi))*cos(x)+10;
    elseif sum(y) == 2
        f = (X(3)-5.1/(4*pi^2)*x^2+5/(pi)*x-6)^2+10*(1-1/(8*pi))*cos(x)+10;
    elseif sum(y) == 3
        f = (X(4)-5.1/(4*pi^2)*x^2+5/(pi)*x-6)^2+10*(1-1/(8*pi))*cos(x)+10;
    end
elseif mode == 7 % Hatman 3
    n_level = 4;
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
    
    if sum(y) == 0
        x(3) = X(1);
        T = 0;
        for i = 1:4
            K = 0;
            for j = 1:3
                K = K + a(j,i)*(x(j)-p(j,i))^2;
            end
            T = T +c(i)*exp(-K);
        end
        f = -T;
    elseif sum(y) == 1
        x(3) = X(2);
        T = 0;
        for i = 1:4
            K = 0;
            for j = 1:3
                K = K + a(j,i)*(x(j)-p(j,i))^2;
            end
            T = T +c(i)*exp(-K);
        end
        f = -T;
    elseif sum(y) == 2
        x(3) = X(3);
        T = 0;
        for i = 1:4
            K = 0;
            for j = 1:3
                K = K + a(j,i)*(x(j)-p(j,i))^2;
            end
            T = T +c(i)*exp(-K);
        end
        f = -T;
    elseif sum(y) == 3
        x(3) = X(4);
        T = 0;
        for i = 1:4
            K = 0;
            for j = 1:3
                K = K + a(j,i)*(x(j)-p(j,i))^2;
            end
            T = T +c(i)*exp(-K);
        end
        f = -T;
    end
elseif mode == 8 % HArtman6
    dimen = 6;
    n_level = 4;
    xl = 0;
    xu = 1;
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
    
    if sum(y) == 0
        x(dimen) = X(1);
        T = 0;
        for i = 1:4
            K = 0;
            for j = 1:6
                K = K + a(j,i)*(x(j)-p(j,i))^2;
            end
            T = T +c(i)*exp(-K);
        end
        f = -T;
    elseif sum(y) == 1
        x(dimen) = X(2);
        T = 0;
        for i = 1:4
            K = 0;
            for j = 1:3
                K = K + a(j,i)*(x(j)-p(j,i))^2;
            end
            T = T +c(i)*exp(-K);
        end
        f = -T;
    elseif sum(y) == 2
        x(dimen) = X(3);
        T = 0;
        for i = 1:4
            K = 0;
            for j = 1:3
                K = K + a(j,i)*(x(j)-p(j,i))^2;
            end
            T = T +c(i)*exp(-K);
        end
        f = -T;
    elseif sum(y) == 3
        x(dimen) = X(4);
        T = 0;
        for i = 1:4
            K = 0;
            for j = 1:3
                K = K + a(j,i)*(x(j)-p(j,i))^2;
            end
            T = T +c(i)*exp(-K);
        end
        f = -T;
    end
elseif mode == 9 %Perm6
    dimen = 6;
    n_level = 4;
    xl = -6;
    xu = 6;
    beta = 60;
    X = zeros(n_level,1);
    for i = 1:n_level
        X(i) = xl + (i-1)*(xu - xl)/(n_level-1);
    end
    
    if sum(y) == 0
        x(dimen) = X(1);
        T = 0;
        for k = 1:6
            S = 0;
            for i = 1:6
                S = S + ((i+1)^k+beta)*(x(i)/(i+1)^k-1);
            end
            T = T + S^2+1000;
        end
        f = T;
    elseif sum(y) == 1
        x(dimen) = X(2);
        T = 0;
        for k = 1:6
            S = 0;
            for i = 1:6
                S = S + ((i+1)^k+beta)*(x(i)/(i+1)^k-1);
            end
            T = T + S^2+1000;
        end
        f = T;
    elseif sum(y) == 2
        x(dimen) = X(3);
        T = 0;
        for k = 1:6
            S = 0;
            for i = 1:6
                S = S + ((i+1)^k+beta)*(x(i)/(i+1)^k-1);
            end
            T = T + S^2+1000;
        end
        f = T;
    elseif sum(y) == 3
        x(dimen) = X(4);
        T = 0;
        for k = 1:6
            S = 0;
            for i = 1:6
                S = S + ((i+1)^k+beta)*(x(i)/(i+1)^k-1);
            end
            T = T + S^2+1000;
        end
        f = T;
    end
elseif mode == 10 % Perm8
    dimen = 8;
    n_level = 4;
    xl = -1;
    xu = 1;
    beta = 100;
    X = zeros(n_level,1);
    for i = 1:n_level
        X(i) = xl + (i-1)*(xu - xl)/(n_level-1);
    end
    
    if sum(y) == 0
        x(dimen) = X(1);
        T = 0;
        for k = 1:8
            S = 0;
            for i = 1:8
                S = S + ((i+1)+beta)*(x(i)^k -(1/(i+1))^k);
            end
            T = T + S^2+1000;
        end
        f = T;
    elseif sum(y) == 1
        x(dimen) = X(2);
        T = 0;
        for k = 1:8
            S = 0;
            for i = 1:8
                S = S + ((i+1)+beta)*(x(i)^k -(1/(i+1))^k);
            end
            T = T + S^2+1000;
        end
        f = T;
    elseif sum(y) == 2
        x(dimen) = X(3);
        T = 0;
        for k = 1:8
            S = 0;
            for i = 1:8
                S = S + ((i+1)+beta)*(x(i)^k -(1/(i+1))^k);
            end
            T = T + S^2+1000;
        end
        f = T;
    elseif sum(y) == 3
        x(dimen) = X(4);
        T = 0;
        for k = 1:8
            S = 0;
            for i = 1:8
                S = S + ((i+1)+beta)*(x(i)^k -(1/(i+1))^k);
            end
            T = T + S^2+1000;
        end
        f = T;
    end
 elseif mode == 11 % Sporttournament
        n_level = 4;
        xl = 0;
        xu = 1;
        X = zeros(n_level,1);
        for i = 1:n_level
            X(i) = xl + (i-1)*(xu - xl)/(n_level-1);
        end
        
        if sum(y) == 0
            f = 2*x(1)*x(3)-2*x(1)+2*x(3)+2*x(1)*x(7)-2*x(7) ...
                +2*x(2)*x(5)-2*x(2)-2*x(5)+2*x(2)*x(10) - ...
                4*x(10)-2*x(3)*x(4)+2*x(4)-2*x(3)*x(12) - ...
                2*x(3)*x(14)-2*x(4)*x(5)+2*x(4)*x(9)-2*x(9) - ...
                2*x(4)*X(1)+2*x(5)*x(6)-2*x(6)+2*x(5)*x(8) - ...
                2*x(8)+2*x(6)*x(9)-2*x(7)*x(8)+2*x(7)*x(12)+ ...
                2*x(7)*x(13)+2*x(8)*x(10)+2*x(8)*X(1) + ...
                2*x(9)*x(11)-2*x(11)-2*x(9)*x(12)+2*x(10)*x(11) + ...
                2*x(10)*x(12)-2*x(13)*X(1)+2*x(14)*X(1);
        elseif sum(y) == 1
            f = 2*x(1)*x(3)-2*x(1)+2*x(3)+2*x(1)*x(7)-2*x(7) ...
                +2*x(2)*x(5)-2*x(2)-2*x(5)+2*x(2)*x(10) - ...
                4*x(10)-2*x(3)*x(4)+2*x(4)-2*x(3)*x(12) - ...
                2*x(3)*x(14)-2*x(4)*x(5)+2*x(4)*x(9)-2*x(9) - ...
                2*x(4)*X(2)+2*x(5)*x(6)-2*x(6)+2*x(5)*x(8) - ...
                2*x(8)+2*x(6)*x(9)-2*x(7)*x(8)+2*x(7)*x(12)+ ...
                2*x(7)*x(13)+2*x(8)*x(10)+2*x(8)*X(2) + ...
                2*x(9)*x(11)-2*x(11)-2*x(9)*x(12)+2*x(10)*x(11) + ...
                2*x(10)*x(12)-2*x(13)*X(2)+2*x(14)*X(2);
        elseif sum(y) == 2
            f = 2*x(1)*x(3)-2*x(1)+2*x(3)+2*x(1)*x(7)-2*x(7) ...
                +2*x(2)*x(5)-2*x(2)-2*x(5)+2*x(2)*x(10) - ...
                4*x(10)-2*x(3)*x(4)+2*x(4)-2*x(3)*x(12) - ...
                2*x(3)*x(14)-2*x(4)*x(5)+2*x(4)*x(9)-2*x(9) - ...
                2*x(4)*X(3)+2*x(5)*x(6)-2*x(6)+2*x(5)*x(8) - ...
                2*x(8)+2*x(6)*x(9)-2*x(7)*x(8)+2*x(7)*x(12)+ ...
                2*x(7)*x(13)+2*x(8)*x(10)+2*x(8)*X(3) + ...
                2*x(9)*x(11)-2*x(11)-2*x(9)*x(12)+2*x(10)*x(11) + ...
                2*x(10)*x(12)-2*x(13)*X(3)+2*x(14)*X(3);
        elseif sum(y) == 3
            f = 2*x(1)*x(3)-2*x(1)+2*x(3)+2*x(1)*x(7)-2*x(7) ...
                +2*x(2)*x(5)-2*x(2)-2*x(5)+2*x(2)*x(10) - ...
                4*x(10)-2*x(3)*x(4)+2*x(4)-2*x(3)*x(12) - ...
                2*x(3)*x(14)-2*x(4)*x(5)+2*x(4)*x(9)-2*x(9) - ...
                2*x(4)*X(4)+2*x(5)*x(6)-2*x(6)+2*x(5)*x(8) - ...
                2*x(8)+2*x(6)*x(9)-2*x(7)*x(8)+2*x(7)*x(12)+ ...
                2*x(7)*x(13)+2*x(8)*x(10)+2*x(8)*X(4) + ...
                2*x(9)*x(11)-2*x(11)-2*x(9)*x(12)+2*x(10)*x(11) + ...
                2*x(10)*x(12)-2*x(13)*X(4)+2*x(14)*X(4);
        end
   
elseif mode == 12 % Extended Hartman 3 - 8 binaries -> 36 levels
    n_level = 36;
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
    for i = 1:n_level
        if myisrotation(y,Neck(i,:))
            x(3) = X(i);
            T = 0;
            for k = 1:4
                K = 0;
                for j = 1:3
                    K = K + a(j,k)*(x(j)-p(j,k))^2;
                end
                T = T +c(k)*exp(-K);
            end
            f = -T;
        end
    end
elseif mode == 13 % Extended Perm 8 - 8 binaries -> 36 levels
    dimen = 8;
    n_level = 36;
    xl = -1;
    xu = 1;
    beta = 100;
    X = zeros(n_level,1);
    for i = 1:n_level
        X(i) = xl + (i-1)*(xu - xl)/(n_level-1);
    end
    for i = 1:n_level
        if myisrotation(y,Neck(i,:))
            x(dimen) = X(i);
            T = 0;
            for k = 1:8
                S = 0;
                for j = 1:8
                    S = S + ((j+1)+beta)*(x(j)^k -(1/(j+1))^k);
                end
                T = T + S^2+1000;
            end
            f = T;
        end
    end
elseif mode == 14 % Extended Branin 8 binaries -> 36 levels
    n_level = 36;
    xl = -10;
    xu = 10;
    X = zeros(n_level,1);
    for i = 1:n_level
        X(i) = xl + (i-1)*(xu - xl)/(n_level-1);
    end
    for i = 1:n_level
        if myisrotation(y,Neck(i,:))
            f = (X(i)-5.1/(4*pi^2)*x^2+5/(pi)*x-6)^2+10*(1-1/(8*pi))*cos(x)+10;
        end
    end
elseif mode == 15 % Extended Branin 10 binaries -< 108 levels
    n_level = 108;
    xl = -10;
    xu = 10;
    X = zeros(n_level,1);
    for i = 1:n_level
        X(i) = xl + (i-1)*(xu - xl)/(n_level-1);
    end
    for i = 1:n_level
        if myisrotation(y,Neck(i,:))
            f = (X(i)-5.1/(4*pi^2)*x^2+5/(pi)*x-6)^2+10*(1-1/(8*pi))*cos(x)+10;
        end
    end
end           
end
