%% NOT FINISH
%***********************************************************%
%              >> IN THE NAME OF GOD <<                     %
% Pushover Analysis of 2nd Linear truss                     %
% subjected to constant axial and and lateral displacement  %
% Checking the analysis by Displacement - Large displecement%
%-----------------------------------------------------------%
%     This program is written by salar delavar ghashghaei   %  
%          E-mail:salar.d.ghashghaei@gmail.com              %
%***********************************************************%
clear all;close all;clc
% Define Parameters in unit: mm,kN
P1=-10; % [kN]
P2=-50; % [kN]
D1=10;% [mm] Initial Displacement  [DOF (1)] Incremantal Displacement
D2=-.01;% [mm] Initial Displacement  [DOF (2)] Incremantal Displacement
D1max=150; % [mm] Maximum displacement [DOF (1)]
D2max=500; % [mm] Maximum displacement [DOF (2)]
XY1i=[250 500]; % [x y] Point 1 Coordinate
XY2i=[0 0]; % [x y] Point 2 Coordinate
XY3i=[260 0]; % [x y] Point 3 Coordinate
XY4i=[500 0]; % [x y] Point 4 Coordinate
A1 = (50)^2; % [mm^2]
A2 = (10)^2; % [mm^2]
A3 = (50)^2; % [mm^2]
E=200;% [N/mm^2] Modulus of elasticity
m = 1000; % number of calculation
itermax = 5000;% maximum number of iterations
tolerance = 1e-9; % specified tolerance for convergence
EA1 = E*A1; % [kN]
EA2 = E*A2; % [kN]
EA3 = E*A3; % [kN]
u = zeros(1,1);% initial guess value
%%% monitor cpu time
starttime = cputime;

%% Large Displacement Analysis [DOF(6)]
disp('########################################');
disp('# Large Displacement Analysis [DOF(1)] #');
disp('########################################');
% Gradually increase the applied displacement
for i=1:m
up=D1*i;
XY1=[XY1i(1)+up XY1i(2)+u(1)]; % [x y] Point 1 Coordinate
XY2=[XY2i(1) XY2i(2)]; % [x y] Point 2 Coordinate
XY3=[XY3i(1) XY3i(2)]; % [x y] Point 3 Coordinate
XY4=[XY4i(1) XY4i(2)]; % [x y] Point 4 Coordinate
L1 = ((XY1(1)-XY2(1))^2+(XY1(2)-XY2(2))^2)^.5; % [mm]
L2 = ((XY1(1)-XY3(1))^2+(XY1(2)-XY3(2))^2)^.5; % [mm]
L3 = ((XY1(1)-XY4(1))^2+(XY1(2)-XY4(2))^2)^.5; % [mm]
lanX1=(XY2(1)-XY1(1))/L1;lanY1=(XY2(2)-XY1(2))/L1;
lanX2=(XY3(1)-XY1(1))/L2;lanY2=(XY3(2)-XY1(2))/L2;
lanX3=(XY4(1)-XY1(1))/L3;lanY3=(XY4(2)-XY1(2))/L3;
G1=EA1/L1;
G2=EA2/L2;
G3=EA2/L3;
% initial assemble global K matrix
%% NOT DONE
        %  1    2
Kp = [G1*lanX1^2+G2*lanX2^2+G3*lanX3^2 G1*lanX1*lanY1+G2*lanX2*lanY2+G3*lanX3*lanY3;
      G1*lanX1*lanY1+G2*lanX2*lanY2+G3*lanX3*lanY3 G1*lanY1^2+G2*lanY2^2+G3*lanY3^2];
Fii=Kp(:,2)*up;
Kini = [G1*lanX1^2+G2*lanX2^2+G3*lanX3^2];
% Define the applied load
Fi = [P1;P2];
F=Fi-Fii;F=[F(1,1)];
        it = 0; % initialize iteration count
        residual = 100; % initialize residual
      while (residual > tolerance)
XY1=[XY1i(1)+up XY1i(2)+u(1)]; % [x y] Point 1 Coordinate
XY2=[XY2i(1) XY2i(2)]; % [x y] Point 2 Coordinate
XY3=[XY3i(1) XY3i(2)]; % [x y] Point 3 Coordinate
XY4=[XY4i(1) XY4i(2)]; % [x y] Point 4 Coordinate
L1 = ((XY1(1)-XY2(1))^2+(XY1(2)-XY2(2))^2)^.5; % [mm]
L2 = ((XY1(1)-XY3(1))^2+(XY1(2)-XY3(2))^2)^.5; % [mm]
L3 = ((XY1(1)-XY4(1))^2+(XY1(2)-XY4(2))^2)^.5; % [mm]
lanX1=(XY2(1)-XY1(1))/L1;lanY1=(XY2(2)-XY1(2))/L1;
lanX2=(XY3(1)-XY1(1))/L2;lanY2=(XY3(2)-XY1(2))/L2;
lanX3=(XY4(1)-XY1(1))/L3;lanY3=(XY4(2)-XY1(2))/L3;
G1=EA1/L1;
G2=EA2/L2;
G3=EA2/L3;
% initial assemble global K matrix
        %  1    2
   K = [G1*lanX1^2+G2*lanX2^2+G3*lanX3^2];        
        f=K*u-F;
        %calculate du1 & du2
        du = Kini^-1 *(-f);
        %Calculate the residual (internal-external force deviation)
        residual = max(abs(du)); % evaluate residual
        it = it + 1; % increment iteration count
        if it == itermax
          fprintf('(-)For increment %1.0f trail iteration reached to Ultimate %1.0f\n',i,it)
             disp('    ## The solution for this step is not converged ##') 
            break
        end
        u = u+du; % update u
      end
              % iteration control
              if it < itermax
              fprintf('(+)Increment %1.0f : It is converged in %1.0f iterations\n',i,it)
              end
lanX1N1(i)=lanX1;lanY1N1(i)=lanY1;
lanX2N1(i)=lanX2;lanY2N1(i)=lanY2;
lanX3N1(i)=lanX3;lanY3N1(i)=lanY3;
%% Force and Dispalcement for each increment
% Internal element force          
         % Displacement Transformation Matrix
        T1 = [lanX1 lanY1 0 0;
              0 0 lanX1 lanX1];
        % Stiffness Matrix for each element
        Kele1 = [G1 -G1;
                -G1 G1]; 
        Fele1 = T1'*Kele1*T1*[0;0;up;u(1)];
        Fele1i1(i) = G1*[-lanX1 -lanY1 lanX1 lanY1]*[0;0;up;u(1)];% Internal Force of element 1
        % Displacement Transformation Matrix
        T2 = [lanX2 lanY2 0 0;
              0 0 lanX2 lanY2];
        Kele2 = [G2 -G2;
                -G2 G2];        
        Fele2 = T2'*Kele2*T2*[0;0;up;u(1)];
        Fele2i1(i) = G2*[-lanX2 -lanY2 lanX2 lanY2]*[0;0;up;u(1)];% Internal Force of element 2
        % Displacement Transformation Matrix
        T3 = [lanX3 lanY3 0 0;
              0 0 lanX3 lanY3];
        Kele3 = [G3 -G3;
                -G3 G3];        
        Fele3 = T3'*Kele3*T3*[0;0;up;u(1)];
        Fele3i1(i) = G3*[-lanX3 -lanY3 lanX3 lanY3]*[0;0;up;u(1)];% Internal Force of element 3
% Force and Dispalcement for each increment
    U1N(i) = up;
    U2N(i) = u(1);
    INT_N_f1(i) = Fele1(1);% Interal force of element 1
    INT_N_f2(i) = Fele2(1);% Interal force of element 2
    INT_N_f3(i) = Fele3(1);% Interal force of element 3
    TBSN1(i)=-[INT_N_f1(i)+INT_N_f2(i)+INT_N_f3(i)]; % Total Base Shear
    if abs(up) >= D1max;disp('      ## Displacement at [DOF (1)] reached to Ultimate Displacement ##');break;end
end
Do1=[U1N'];Fo1=[TBSN1'];
disp('= Large Displacement Analysis [DOF(1)] =');
disp(' X-displacement(D1)   Y-displacement(D2)');
disp('----------------------------------------')
disp([U1N' U2N'])
disp('+=================================================+');
disp('  Axial Load(f1)   Axial Load(f2))   Axial Load(f3)');
disp('---------------------------------------------------');
disp([INT_N_f1' INT_N_f2' INT_N_f3']);
disp('+======================================+');
disp('  Y-displacement(D1)   Base Reaction ');
disp('----------------------------------------');
disp([U1N' TBSN1']);
%% SAP2000 Analysis
sapx=[0
20.005083
40.833751
61.690781
82.512955
103.299041
124.052471
144.776603
165.474413
186.148534
206.801308];
sapy=-[0
730.292
1491.402
2254.402
3016.945
3778.98
4540.63
5302.013
6063.234
6824.387
7585.556];
%% Large Displacement Analysis [DOF(5)]
disp('########################################');
disp('# Large Displacement Analysis [DOF(2)] #');
disp('########################################');
u = zeros(1,1);% initial guess value
% Gradually increase the applied displacement
for i=1:m
up=D2*i;
XY1=[XY1i(1)+u(1) XY1i(2)+up]; % [x y] Point 1 Coordinate
XY2=[XY2i(1) XY2i(2)]; % [x y] Point 2 Coordinate
XY3=[XY3i(1) XY3i(2)]; % [x y] Point 3 Coordinate
XY4=[XY4i(1) XY4i(2)]; % [x y] Point 4 Coordinate
L1 = ((XY1(1)-XY2(1))^2+(XY1(2)-XY2(2))^2)^.5; % [mm]
L2 = ((XY1(1)-XY3(1))^2+(XY1(2)-XY3(2))^2)^.5; % [mm]
L3 = ((XY1(1)-XY4(1))^2+(XY1(2)-XY4(2))^2)^.5; % [mm]
lanX1=(XY2(1)-XY1(1))/L1;lanY1=(XY2(2)-XY1(2))/L1;
lanX2=(XY3(1)-XY1(1))/L2;lanY2=(XY3(2)-XY1(2))/L2;
lanX3=(XY4(1)-XY1(1))/L3;lanY3=(XY4(2)-XY1(2))/L3;
G1=EA1/L1;
G2=EA2/L2;
G3=EA2/L3;
% initial assemble global K matrix
        %  1    2
Kp = [G1*lanX1^2+G2*lanX2^2+G3*lanX3^2 G1*lanX1*lanY1+G2*lanX2*lanY2+G3*lanX3*lanY3;
      G1*lanX1*lanY1+G2*lanX2*lanY2+G3*lanX3*lanY3 G1*lanY1^2+G2*lanY2^2+G3*lanY3^2];
Fii=Kp(:,1)*up;
Kini = [G1*lanY1^2+G2*lanY2^2+G3*lanY3^2];
% Define the applied load
Fi = [P1;P2];
F=Fi-Fii;F=[F(2,1)];
        it = 0; % initialize iteration count
        residual = 100; % initialize residual
      while (residual > tolerance)
XY1=[XY1i(1)+u(1) XY1i(2)+up]; % [x y] Point 1 Coordinate
XY2=[XY2i(1) XY2i(2)]; % [x y] Point 2 Coordinate
XY3=[XY3i(1) XY3i(2)]; % [x y] Point 3 Coordinate
XY4=[XY4i(1) XY4i(2)]; % [x y] Point 4 Coordinate
L1 = ((XY1(1)-XY2(1))^2+(XY1(2)-XY2(2))^2)^.5; % [mm]
L2 = ((XY1(1)-XY3(1))^2+(XY1(2)-XY3(2))^2)^.5; % [mm]
L3 = ((XY1(1)-XY4(1))^2+(XY1(2)-XY4(2))^2)^.5; % [mm]
lanX1=(XY2(1)-XY1(1))/L1;lanY1=(XY2(2)-XY1(2))/L1;
lanX2=(XY3(1)-XY1(1))/L2;lanY2=(XY3(2)-XY1(2))/L2;
lanX3=(XY4(1)-XY1(1))/L3;lanY3=(XY4(2)-XY1(2))/L3;
G1=EA1/L1;
G2=EA2/L2;
G3=EA2/L3;
% initial assemble global K matrix
        %  1    2
   K = [G1*lanY1^2+G2*lanY2^2+G3*lanY3^2];        
        f=K*u-F;
        %calculate du1 & du2
        du = Kini^-1 *(-f);
        %Calculate the residual (internal-external force deviation)
        residual = max(abs(du)); % evaluate residual
        it = it + 1; % increment iteration count
        if it == itermax
          fprintf('(-)For increment %1.0f trail iteration reached to Ultimate %1.0f\n',i,it)
             disp('    ## The solution for this step is not converged ##') 
            break
        end
        u = u+du; % update u
      end
              % iteration control
              if it < itermax
              fprintf('(+)Increment %1.0f : It is converged in %1.0f iterations\n',i,it)
              end
lanX1L1(i)=lanX1;lanY1L1(i)=lanY1;
lanX2L1(i)=lanX2;lanY2L1(i)=lanY2;
lanX3L1(i)=lanX3;lanY3L1(i)=lanY3;
%% Force and Dispalcement for each increment
% Internal element force          
         % Displacement Transformation Matrix
        T1 = [lanX1 lanY1 0 0;
              0 0 lanX1 lanY1];
        % Stiffness Matrix for each element
        Kele1 = [G1 -G1;
                -G1 G1]; 
        Fele1 = T1'*Kele1*T1*[0;0;u(1);up];
        Fele1i2(i) = G1*[-lanX1 -lanY1 lanX1 lanY1]*[0;0;u(1);up];% Internal Force of element 1
        % Displacement Transformation Matrix
        T2 = [lanX2 lanY2 0 0;
              0 0 lanX2 lanY2];
        Kele2 = [G2 -G2;
                -G2 G2];        
        Fele2 = T2'*Kele2*T2*[0;0;u(1);up];
        Fele2i2(i) = G2*[-lanX2 -lanY2 lanX2 lanY2]*[0;0;u(1);up];% Internal Force of element 2
        % Displacement Transformation Matrix
        T3 = [lanX3 lanY3 0 0;
              0 0 lanX3 lanY3];
        Kele3 = [G3 -G3;
                -G3 G3];        
        Fele3 = T3'*Kele3*T3*[0;0;u(1);up];
        Fele3i2(i) = G3*[-lanX3 -lanY3 lanX3 lanY3]*[0;0;u(1);up];% Internal Force of element 1        
% Force and Dispalcement for each increment
    U1L(i) = u(1);
    U2L(i) = up;
    INT_L_f1(i) = Fele1(1);% Interal force of element 1
    INT_L_f2(i) = Fele2(1);% Interal force of element 2
    INT_L_f3(i) = Fele3(1);% Interal force of element 3
    TBSN2(i)=-[INT_L_f1(i)+INT_L_f2(i)+INT_L_f3(i)]; % Total Base Shear
end
Do2=[U2L'];Fo2=[TBSN2'];
disp('= Large Displacement Analysis [DOF(2)] =');
disp(' X-displacement(D1)   Y-displacement(D2)');
disp('----------------------------------------')
disp([U1L' U2L'])
disp('+=================================================+');
disp('  Axial Load(f1)   Axial Load(f2))   Axial Load(f3)');
disp('---------------------------------------------------');
disp([INT_L_f1' INT_L_f2' INT_L_f3']);
disp('+======================================+');
disp('  Y-displacement(D2)   Base Reaction ');
disp('----------------------------------------');
disp([U1L' TBSN2']);
disp('+======================================+');
%% imaging
figure(1)
IMAGE=imread('PushoverLateralTruss02.jpg');
image(IMAGE);axis image;axis off;
figure(2)
p1=plot(lanX1N1,lanY1N1,'-.black',lanX2N1,lanY2N1,'-.blue');grid on;set(p1,'LineWidth',3);
legend('Angle of element 1','Angle of element 2','Location','NorthEastOutside');
xlabel('Angle of element during the incremental Loading [cos (a)] (rad)');ylabel('Angle of element during the incremental Loading [sin (a)] (rad)]');
title('[DOF(1)]: Angle of elements during the incremental displacement','color','b')
figure(3)
p1=plot(U1N,U2N,'oblack');grid on;set(p1,'LineWidth',3);
xlabel('Displacement-X (mm) [DOF(1)]');ylabel('Displacement-Y (mm) [DOF(2)]');
title('[DOF(1)]: Large Displacement Theory - Displacement-X and Y during the incremental displacement[DOF(6)]','color','b');
figure(4)
p1=plot(Do1,Fo1,'r--');grid on;set(p1,'LineWidth',3);
xlabel('Displacement (mm) [DOF(1)]');ylabel('Base Reaction (kN) [DOF(4)+[DOF(6)+[DOF(8)]');
title('Base Reaction-Displacement of truss during the incremental displacement[DOF(1)]','color','b');
legend('Analysis','SAP2000','Location','NorthEastOutside');
figure(5)
p1=plot(Do1,Fele1i1,Do1,Fele2i1,'--',Do1,Fele3i1);grid on;set(p1,'LineWidth',3);
xlabel('Displacement (mm) [DOF(1)]');ylabel('Axial Load of element (kN)');
title('Element Axial Load-Displacement of truss during the incremental displacement [DOF(1)] (Positive=Tension) (Negative=Compression)','color','b');
legend('ele. 1','ele. 2','ele. 3','Location','NorthEastOutside');
figure(6)
p1=plot(lanX1L1,lanY1L1,'-.black',lanX2L1,lanY2L1,'-.blue');grid on;set(p1,'LineWidth',3);
legend('Angle of element 1','Angle of element 2','Location','NorthEastOutside');
xlabel('Angle of element during the incremental Loading [cos (a)] (rad)');
ylabel('Angle of element during the incremental Loading [sin (a)] (rad)]');
title('[DOF(2)]: Angle of elements during the incremental displacement','color','b')
figure(7)
p1=plot(U1L,U2L,'oblack');grid on;set(p1,'LineWidth',3);
xlabel('Displacement-X (mm) [DOF(1)]');ylabel('Displacement-Y (mm) [DOF(2)]');
title('[DOF(2)]: Large Displacement Theory - Displacement-X and Y during the incremental displacement','color','b');
figure(8)
p1=plot(Do2,Fo2);grid on;set(p1,'LineWidth',3);
xlabel('Displacement (mm) [DOF(2)]');ylabel('Base Reaction (kN) [DOF(4)+[DOF(6)+[DOF(8)]');
title('Base Reaction-Displacement of truss during the incremental displacement[DOF(2)]','color','b');
figure(9)
p1=plot(Do2,Fele1i2,Do2,Fele2i2,'--',Do2,Fele3i2);grid on;set(p1,'LineWidth',3);
xlabel('Displacement (mm) [DOF(1)]');ylabel('Axial Load of element (kN)');
title('Element Axial Load-Displacement of truss during the incremental displacement [DOF(1)]  (Positive=Tension) (Negative=Compression)','color','b');
legend('ele. 1','ele. 2','ele. 3','Location','NorthEastOutside');