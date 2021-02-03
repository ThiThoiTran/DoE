
% function of testing, cite Luksan & Vlcek (2000)
function [f,g] = test(mode,x,y)
    %% minimax unconstrained funcitons
    %%% NC = 2. NB = 2 (nlevel = 3)
    if mode == 1 % CB2 (2,2)
        if sum(y) == 0
            f = x(1)^2 + x(2)^4;
            g = [2*x(1); 2*x(2)];
        elseif sum(y) == 1
            f = (2 - x(1))^2 + (2 - x(2))^2;
            g = [-2*(2 - x(1)); -2*(2 - x(2))];
        elseif sum(y) == 2
            f = 2*exp(x(2) - x(1));
            g = [-2*exp(x(2) - x(1)); 2*exp(x(2) - x(1))];
        end
   
        %% general nonsmooth functions (unconstrained)
    elseif mode == 2 % CB3 (2,2)
        if sum(y) == 0
            f = x(1)^4 + x(2)^2;
            g = [2*x(1); 2*x(2)];
        elseif sum(y) == 1
            f = (2 - x(1))^2 + (2 - x(2))^2;
            g = [-2*(2 - x(1)); -2*(2 - x(2))];
        elseif sum(y) == 2
            f = 2*exp(x(2) - x(1));
            g = [-2*exp(x(2) - x(1)); 2*exp(x(2) - x(1))];
        end

        %% general nonsmooth functions (unconstrained)
  
    elseif mode == 3 % QL (-1,5)
        if sum(y) == 0
            f = x(1)^2 + x(2)^2;
            g = [2*x(1); 2*x(2)];
        elseif sum(y) == 1
            f = x(1)^2 + x(2)^2 + 10*(-4*x(1) - x(2) + 4);
            g = [];
        elseif sum(y) == 2
            f = x(1)^2 + x(2)^2 + 10*(-x(1) - 2*x(2) + 6);
            g = [];
        end
   elseif mode == 4% WF (3,1)
        if sum(y) == 1
            f = 1/2*(x(1) + (10*x(1))/(x(1) + 0.1) + 2*x(2)^2);
            g = [];
        elseif sum(y) == 0
            f = 1/2*(-x(1) + (10*x(1))/(x(1) + 0.1) + 2*x(2)^2);
            g = [];
        elseif sum(y) == 2
            f = 1/2*(x(1) - (10*x(1))/(x(1) + 0.1) + 2*x(2)^2);
            g = [];
        end
    elseif mode == 5 % EVD52
       if isequal(y,[0 0 0 0]) 
           f = x(1)^2 + x(2)^2 + (x(3) -2)^2; %f = x(1)^2 + x(2)^2 + x(3)^2 -1;
       elseif isequal(y,[1 0 0 0]) || isequal(y,[0 1 0 0])|| isequal(y,[0 0 1 0])...
               || isequal(y,[0 0 0 1]) 
           f = x(1)^2 + x(2)^2 + x(3)^2 -1; %f = x(1)^2 + x(2)^2 + (x(3) -2)^2;
       elseif isequal(y,[1 0 1 0])|| isequal(y,[0 1 0 1])
           f = x(1) + x(2) - x(3) + 1; %f = x(1) + x(2) + x(3) - 1;
       elseif isequal(y,[1 1 0 0])||isequal(y,[0 1 1 0])||isequal(y,[0 0 1 1])...
               ||isequal(y,[1 0 0 1]) 
          f = x(1) + x(2) + x(3) - 1; 
       elseif isequal(y,[1 1 1 0])||isequal(y,[0 1 1 1])||isequal(y,[1 0 1 1])...
               ||isequal(y,[1 1 0 1]) 
            f = x(1)^2 - 9*x(3);  %f = 2*x(1)^3 + 6*x(2)^2 + 2*(5*x(3) - x(1) + 1)^2;
       elseif isequal(y,[1 1 1 1])
           f = 2*x(1)^3 + 6*x(2)^2 + 2*(5*x(3) - x(1) + 1)^2; %f = x(1)^2 - 9*x(3);
       end
       g = [];
   elseif mode == 6 % Pentagon
       if sum(y) == 0
           f = -sqrt((x(1) - x(3))^2 + (x(2) - x(4))^2);
       elseif sum(y) == 1
           f = -sqrt((x(3) - x(5))^2 + (x(4) - x(6))^2);
       elseif sum(y) == 2
           f = -sqrt((x(5) - x(1))^2 + (x(6) - x(2))^2);
       end
       g = [];
        %% NC = 4, NB = 3 (nlevel = 4)
    elseif mode == 7 % Rosen-Suzuki (0,0,0,0)
        if sum(y) == 0
            f = x(1)^2 + x(2)^2 + 2*x(3)^2 + x(4)^2 - 5*x(1) - 5*x(2) - 21*x(3) + 7*x(4);
            g = [];
        elseif sum(y) == 1
            f = x(1)^2 + x(2)^2 + 2*x(3)^2 + x(4)^2 - 5*x(1) - 5*x(2) - 21*x(3) + 7*x(4)...
                + 10*(x(1)^2 + x(2)^2 + x(3)^2 + x(4)^2 + x(1) -x(2) + x(3) -x(4) -8);
            g = [];
        elseif sum(y) == 2
            f = x(1)^2 + x(2)^2 + 2*x(3)^2 + x(4)^2 - 5*x(1) - 5*x(2) - 21*x(3) + 7*x(4)...
                + 10*(x(1)^2 + 2*x(2)^2 + x(3)^2 + 2*x(4)^2 - x(1) -x(4)- 10);
            g = [];
        elseif sum(y) == 3
            f = x(1)^2 + x(2)^2 + 2*x(3)^2 + x(4)^2 - 5*x(1) - 5*x(2) - 21*x(3) + 7*x(4)...
                + 10*(2*x(1)^2 + x(2)^2 + x(3)^2 + 2*x(1) -x(2) -x(4)- 5);
            g = [];
        end
     elseif mode == 8 % DEMBO5
       ref = x(1) + x(2) + x(3);
       a = 833.33252;
       b = -83333.333;
       c = 1250000;
       if sum(y) == 0
            f = ref;
       elseif sum(y) == 1
           f = ref + 10^5*(a*(x(1)^(-1))*x(4)*(x(6)^(-1)) + 100*(x(6)^(-1)) + b*(x(1)^(-1))*(x(6)^(-1)) - 1);
       elseif sum(y) == 2
           f = ref + 10^5*(x(4)*(x(7)^(-1)) + 1250*(x(5)-x(4))*(x(2)^(-1))*(x(7)^(-1))- 1);
       elseif sum(y) == 3
           f = ref + 10^5*(c*(x(3)^(-1))*(x(8)^(-1)) + x(5)*(x(8)^(-1)) - 2500*(x(3)^(-1))*x(5)*(x(8)^(-1))- 1);
       end
       g = [];
       elseif mode == 9 % wong2 (10,4)
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
        g = [];
   elseif mode == 10 % Wong3 (20,6)
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
        g = [];
            
   elseif mode == 11 % MAD1
        if sum(y) == 0
            f = x(1)^2+x(2)^2 + x(1)*x(2) -1;
        elseif sum(y) == 1
            f = sin(x(1));
        elseif sum(y) == 2
            f = -cos(x(2));
        end
        g = [];
    elseif mode == 12 % MAD 2
        if sum(y) == 0
            f = -exp(x(1)-x(2));
        elseif sum(y) == 1
            f = real((exp(x(1) - 1)-exp(-(x(1) - 1)))/2 -1);
        elseif sum(y) == 2
            f = -log(x(2)) - 1;
        end
        g = [];
        
        
    elseif mode == 13 % HS78, take the last variable as Luizzi
        if sum(y) == 0
            f_1 = x(1)^2+x(2)^2+x(3)^2+x(4)^2+X(1)^2-10;        
            f_2 = x(2)*x(3)-5*x(4)*X(1);
            f_3 = x(1)^3 + x(2)^3 + 1;
            f = x(1)*x(2)*x(3)*x(4)*X(1) + 10*(abs(f_1)+abs(f_2)+abs(f_3));
        elseif sum(y) == 1
            f_1 = x(1)^2+x(2)^2+x(3)^2+x(4)^2+X(2)^2-10;        
            f_2 = x(2)*x(3)-5*x(4)*X(2);
            f_3 = x(1)^3 + x(2)^3 + 1;
            f = x(1)*x(2)*x(3)*x(4)*X(2) + 10*(abs(f_1)+abs(f_2)+abs(f_3));
        elseif sum(y) == 2
            f_1 = x(1)^2+x(2)^2+x(3)^2+x(4)^2+X(3)^2-10;        
            f_2 = x(2)*x(3)-5*x(4)*X(3);
            f_3 = x(1)^3 + x(2)^3 + 1;
            f = x(1)*x(2)*x(3)*x(4)*X(3) + 10*(abs(f_1)+abs(f_2)+abs(f_3));
        elseif sum(y) == 3
            f_1 = x(1)^2+x(2)^2+x(3)^2+x(4)^2+X(4)^2-10;        
            f_2 = x(2)*x(3)-5*x(4)*X(4);
            f_3 = x(1)^3 + x(2)^3 + 1;
            f = x(1)*x(2)*x(3)*x(4)*X(4) + 10*(abs(f_1)+abs(f_2)+abs(f_3));
        end
        g = [];
        %% Shift to Hok and schizkovski
    elseif mode == 14;% HS2 (Rosen brock)
        n_level = 4;
        xl = -5;
        xu = 5;
        X = zeros(n_level,1);
        for i = 1:n_level
            X(i) = xl + (i-1)*(xu - xl)/10;
        end
        if sum(y) == 0
            f = 100*(X(1)-x(1)^2)^2 + (1-x(1))^2;
        elseif sum(y) == 1
            f = 100*(X(2)-x(1)^2)^2 + (1-x(1))^2;
        elseif sum(y) == 2
           f = 100*(X(3)-x(1)^2)^2 + (1-x(1))^2;
       elseif sum(y) == 3
           f = 100*(X(4)-x(1)^2)^2 + (1-x(1))^2;
        end
        g = [];
    elseif mode == 15 %HS3
        n_level = 4;
        xl = -5;
        xu = 5;
        X = zeros(n_level,1);
        for i = 1:n_level
            X(i) = xl + (i-1)*(xu - xl)/10;
        end
        if sum(y) == 0
            f = X(1) + 10^(-5)*(X(1) - x(1))^2;
        elseif sum(y) == 1
            f = X(2) + 10^(-5)*(X(2) - x(1))^2;
        elseif sum(y) == 2
            f = X(3) + 10^(-5)*(X(3) - x(1))^2;
       elseif sum(y) == 3
            f = X(4) + 10^(-5)*(X(4) - x(1))^2;
        end
        g = [];
    elseif mode == 16 % HS29log
        n_level = 4;
        xl = -5;
        xu = 5;
        X = zeros(n_level,1);
        for i = 1:n_level
            X(i) = xl + (i-1)*(xu - xl)/10;
        end
        if sum(y) == 0
            f = log10(100*(X(1)-x(1)^2)^2 + (1-x(1))^2);
        elseif sum(y) == 1
            f = log10(100*(X(2)-x(1)^2)^2 + (1-x(1))^2);
        elseif sum(y) == 2
            f = log10(100*(X(3)-x(1)^2)^2 + (1-x(1))^2);
        elseif sum(y) == 3
            f = log10(100*(X(4)-x(1)^2)^2 + (1-x(1))^2);
        end
        g = [];  
    elseif mode == 17 % Branin
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
        g = [];
    elseif mode == 18 %Sixhum
        n_level = 4;
        xl = -2;
        xu = 2;
        X = zeros(n_level,1);
        for i = 1:n_level
            X(i) = xl + (i-1)*(xu - xl)/(n_level-1);
        end
        if sum(y) == 0
            f = (4-2.1*x(1)^2+x(1)^(4))*x(1)^2+x(1)*X(1)+(-4+4*X(1)^2)*X(1)^2;
        elseif sum(y) == 1
            f = (4-2.1*x(1)^2+x(1)^(4))*x(1)^2+x(1)*X(2)+(-4+4*X(2)^2)*X(2)^2;
        elseif sum(y) == 2
            f = (4-2.1*x(1)^2+x(1)^(4))*x(1)^2+x(1)*X(3)+(-4+4*X(3)^2)*X(3)^2;
        elseif sum(y) == 3
            f = (4-2.1*x(1)^2+x(1)^(4))*x(1)^2+x(1)*X(4)+(-4+4*X(4)^2)*X(4)^2;
        end
        g = [];
    elseif mode == 19 % goldstein-price
        n_level = 4;
        xl = -2;
        xu = 2;
        X = zeros(n_level,1);
        for i = 1:n_level
            X(i) = xl + (i-1)*(xu - xl)/(n_level-1);
        end
        if sum(y) == 0
            f = (1+(x(1)+X(1)+1)^2*(19-14*x(1)+3*x(1)^2-14*X(1)+6*x(1)*X(1)+3*X(1)^2))* ...
                (30+(2*x(1)-3*X(1))^2*(18-32*x(1)+12*x(1)^2+48*X(1)-36*x(1)*X(1)+27*X(1)^2));
        elseif sum(y) == 1
            f = (1+(x(1)+X(2)+1)^2*(19-14*x(1)+3*x(1)^2-14*X(2)+6*x(1)*X(2)+3*X(2)^2))* ...
                (30+(2*x(1)-3*X(2))^2*(18-32*x(1)+12*x(1)^2+48*X(2)-36*x(1)*X(2)+27*X(2)^2));
        elseif sum(y) == 2
            f = (1+(x(1)+X(3)+1)^2*(19-14*x(1)+3*x(1)^2-14*X(3)+6*x(1)*X(3)+3*X(3)^2))* ...
                (30+(2*x(1)-3*X(3))^2*(18-32*x(1)+12*x(1)^2+48*X(3)-36*x(1)*X(3)+27*X(3)^2));
        elseif sum(y) == 3
            f = (1+(x(1)+X(4)+1)^2*(19-14*x(1)+3*x(1)^2-14*X(4)+6*x(1)*X(4)+3*X(4)^2))* ...
                (30+(2*x(1)-3*X(4))^2*(18-32*x(1)+12*x(1)^2+48*X(4)-36*x(1)*X(4)+27*X(4)^2));
        end
        g = [];
    elseif mode == 20 % Hartman 3
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
        g = [];
    elseif mode == 21 % hartman 6
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
        g = [];
    elseif mode == 22 %Shekel 7
        dimen = 4;
        n_level = 4;
        xl = 0;
        xu = 10;
        X = zeros(n_level,1);
        for i = 1:n_level
            X(i) = xl + (i-1)*(xu - xl)/(n_level-1);
        end
        
        a = [ [4.0, 1.0, 8.0, 6.0, 3.0, 2.0, 5.0];
            [4.0, 1.0, 8.0, 6.0, 7.0, 9.0, 5.0];
            [4.0, 1.0, 8.0, 6.0, 3.0, 2.0, 3.0];
            [4.0, 1.0, 8.0, 6.0, 7.0, 9.0, 3.0] ];
        c = [0.1, 0.2, 0.2, 0.4, 0.4, 0.6, 0.3];
        if sum(y) == 0
            x(dimen) = X(1);
            T = 0;
            for j = 1:7
                S = 0;
                for i = 1:4
                    S = S + (x(i)-a(i,j))^2;
                end
                T = T + 1/(S+c(j));
            end
            f = -T;
        elseif sum(y) == 1
            x(dimen) = X(2);
            T = 0;
            for j = 1:7
                S = 0;
                for i = 1:4
                    S = S + (x(i)-a(i,j))^2;
                end
                T = T + 1/(S+c(j));
            end
            f = -T;
        elseif sum(y) == 2
            x(dimen) = X(3);
            T = 0;
            for j = 1:7
                S = 0;
                for i = 1:4
                    S = S + (x(i)-a(i,j))^2;
                end
                T = T + 1/(S+c(j));
            end
            f = -T;
        elseif sum(y) == 3
            x(dimen) = X(4);
            T = 0;
            for j = 1:7
                S = 0;
                for i = 1:4
                    S = S + (x(i)-a(i,j))^2;
                end
                T = T + 1/(S+c(j));
            end
            f = -T;
        end
         g = [];
    elseif mode == 23 %Shekel 10
        dimen = 4;
        n_level = 4;
        xl = 0;
        xu = 10;
        X = zeros(n_level,1);
        for i = 1:n_level
            X(i) = xl + (i-1)*(xu - xl)/(n_level-1);
        end
        
        a = [ [4.0, 1.0, 8.0, 6.0, 3.0, 2.0, 5.0, 8.0, 6.0, 7.0];
            [4.0, 1.0, 8.0, 6.0, 7.0, 9.0, 5.0, 1.0, 2.0, 3.6];
            [4.0, 1.0, 8.0, 6.0, 3.0, 2.0, 3.0, 8.0, 6.0, 7.0];
            [4.0, 1.0, 8.0, 6.0, 7.0, 9.0, 3.0, 1.0, 2.0, 3.6] ];
        c = [0.1, 0.2, 0.2, 0.4, 0.4, 0.6, 0.3, 0.7, 0.5, 0.5];
        if sum(y) == 0
            x(dimen) = X(1);
            T = 0;
            for j = 1:10
                S = 0;
                for i = 1:4
                    S = S + (x(i)-a(i,j))^2;
                end
                T = T + 1/(S+c(j));
            end
            f = -T;
        elseif sum(y) == 1
            x(dimen) = X(2);
            T = 0;
            for j = 1:10
                S = 0;
                for i = 1:4
                    S = S + (x(i)-a(i,j))^2;
                end
                T = T + 1/(S+c(j));
            end
            f = -T;
        elseif sum(y) == 2
            x(dimen) = X(3);
            T = 0;
            for j = 1:10
                S = 0;
                for i = 1:4
                    S = S + (x(i)-a(i,j))^2;
                end
                T = T + 1/(S+c(j));
            end
            f = -T;
        elseif sum(y) == 3
            x(dimen) = X(4);
            T = 0;
            for j = 1:10
                S = 0;
                for i = 1:4
                    S = S + (x(i)-a(i,j))^2;
                end
                T = T + 1/(S+c(j));
            end
            f = -T;
        end
        g = [];
    elseif mode == 24 % ex8_11
        n_level = 4;
        xl = -1;
        xu = 1;
        X = zeros(n_level,1);
        for i = 1:n_level
            X(i) = xl + (i-1)*(xu - xl)/(n_level-1);
        end
        
        if sum(y) == 0
            f = cos(x)*sin(X(1)) - x/(X(1)^2+1);
        elseif sum(y) == 1
            f = cos(x)*sin(X(2)) - x/(X(2)^2+1);
        elseif sum(y) == 2
            f = cos(x)*sin(X(3)) - x/(X(3)^2+1);
        elseif sum(y) == 3
            f = cos(x)*sin(X(4)) - x/(X(4)^2+1);
        end
       g = [];
    elseif mode == 25 % ex8_14
        n_level = 4;
        xl = -5;
        xu = 2;
        X = zeros(n_level,1);
        for i = 1:n_level
            X(i) = xl + (i-1)*(xu - xl)/(n_level-1);
        end
        
        if sum(y) == 0
            f = 12*x^2-6.3*x^4+x^6-6*x*X(1)+6*X(1)^2;
        elseif sum(y) == 1
            f = 12*x^2-6.3*x^4+x^6-6*x*X(2)+6*X(2)^2;
        elseif sum(y) == 2
            f = 12*x^2-6.3*x^4+x^6-6*x*X(3)+6*X(3)^2;
        elseif sum(y) == 3
            f = 12*x^2-6.3*x^4+x^6-6*x*X(4)+6*X(4)^2;
        end
        g = [];
    elseif mode == 26 % Perm6
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
        
        g = [];
    elseif mode == 27 % perm 8
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
        
        g = [];
    elseif mode == 28 % Sporttournament
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
        g = [];
        %% Toy problem
    elseif mode == 29
        Nblades = 6;
        m = 0.0114; % unit mass
        kb = 430.3; % stifness blades
        kc = 45.3;% blades-blades stifness
        %kc = 8.606;
        F0 = 1; % initial force
        d = 0.143; % damping value
        OE = 6; % order of exitation
        dev_area = 0.01;%derivation
        dev= 0.1;
        mis = 0.1; %mistuning kb
        %define the shifted phase
        alpha = 2*pi/Nblades;
        % define the left handsid matrix T A = F
        T = zeros(Nblades, Nblades);
        for i = 1: Nblades
            T(i,i)= -m*x^2+(1 + y(i)*(-mis) + (1-y(i))*mis )*kb + 1i*x*d+2*kc;
            % T(i,i)= -m*x^2+(1 + y(i)*mis )*kb + 1i*x*d+2*kc;
        end
        for i = 1:Nblades-1
            T(i,i+1) = -kc;
            T(i+1,i) = -kc;
        end

        T(1,end) = -kc;
        T(end,1) = -kc;
        %Define the right handside function F
        F = zeros(Nblades,1);
        for i = 1:Nblades
            F(i) = F0*exp(1i*(i-1)*alpha*OE);
        end
        % Define cost function
        Ampli = T^(-1)*F;
        f = norm(Ampli, inf) ;
        g = [];
    end
end
        
    
   
