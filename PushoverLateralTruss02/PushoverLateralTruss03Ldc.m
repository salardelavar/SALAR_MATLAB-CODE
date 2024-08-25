%***********************************************************%
%              >> IN THE NAME OF GOD <<                     %
% Pushover Analysis of 1st Linear truss                     %
% subjected to constant axial and and lateral displacement  %
% Checking the analysis by Displacement - Large displacement%
% ABAQUS anad SEISMSOSTRUCT verification                    %
%-----------------------------------------------------------%
%     This program is written by salar delavar ghashghaei   %  
%          E-mail:salar.d.ghashghaei@gmail.com              %
%***********************************************************%
clear all;close all;clc
% Define Parameters in unit: mm,kN
P5=0; % [kN]
P6=0; % [kN]
D5=1;% [mm] Initial Displacement  [DOF (5)] Incremantal Displacement
D6=1;% [mm] Initial Displacement  [DOF (6)] Incremantal Displacement
D5max=100; % [mm] Maximum displacement [DOF (5)]
D6max=500; % [mm] Maximum displacement [DOF (6)]
XY1i=[0 0]; % [x y] Point 1 Coordinate
XY2i=[500 0]; % [x y] Point 2 Coordinate
XY3i=[500 250]; % [x y] Point 3 Coordinate
A1 = 3.1415*(50)^2/4; % [mm^2]
A2 = 3.1415*(1000)^2/4; % [mm^2]
E=200;% [N/mm^2] Modulus of elasticity
m = 100; % number of calculation
itermax = 5000;% maximum number of iterations
tolerance = 1e-12; % specified tolerance for convergence
L1i=(((XY3i(1)-XY1i(1))^2+(XY3i(2)-XY1i(2))^2)^.5);
L2i=(((XY3i(1)-XY2i(1))^2+(XY3i(2)-XY2i(2))^2)^.5);
EA1 = E*A1; % [kN]
EA2 = E*A2; % [kN]
u = zeros(1,1);% initial guess value
%%% monitor cpu time
starttime = cputime;

%% Large Displacement Analysis [DOF(6)]
disp('########################################');
disp('# Large Displacement Analysis [DOF(6)] #');
disp('########################################');
% Gradually increase the applied displacement
for i=1:m
up=D6*i; 
XY1=[XY1i(1) XY1i(2)]; % [x y] Point 1 Coordinate
XY2=[XY2i(1) XY2i(2)]; % [x y] Point 2 Coordinate
XY3=[XY3i(1)+u(1) XY3i(2)+up]; % [x y] Point 3 Coordinate
L1 = ((XY3(1)-XY1(1))^2+(XY3(2)-XY1(2))^2)^.5; % [mm]
L2 = ((XY3i(1)-XY2i(1))^2+(XY3i(2)-XY2i(2))^2)^.5; % [mm]       
lanX1=(XY3(1)-XY1(1))/L1;lanY1=(XY3(2)-XY3(2))/L1;
lanX2=(XY3(1)-XY2(1))/L2;lanY2=(XY3(2)-XY2(2))/L2;    
G1=EA1/L1;
G2=EA2/L2;
% initial assemble global K matrix
        %  5    6
Kp = [G1*lanX1^2+G2*lanX2^2 G1*lanX1*lanY1+G2*lanX2*lanY2;
      G1*lanX1*lanY1+G2*lanX2*lanY2 G1*lanY1^2+G2*lanY2^2];
Fii=Kp(:,2)*up;
Kini = [G1*lanX1^2+G2*lanX2^2];
% Define the applied load
Fi = [P5;P6];
F=Fi-Fii;F=[F(1,1)];
        it = 0; % initialize iteration count
        residual = 100; % initialize residual
while (residual > tolerance)
XY1=[XY1i(1) XY1i(2)]; % [x y] Point 1 Coordinate
XY2=[XY2i(1) XY2i(2)]; % [x y] Point 2 Coordinate
XY3=[XY3i(1)+u(1) XY3i(2)+up]; % [x y] Point 3 Coordinate          
L1 = ((XY3(1)-XY1(1))^2+(XY3(2)-XY1(2))^2)^.5; % [mm]
L2 = ((XY3i(1)-XY2i(1))^2+(XY3i(2)-XY2i(2))^2)^.5; % [mm]       
lanX1=(XY3(1)-XY1(1))/L1;lanY1=(XY3(2)-XY3(2))/L1;
lanX2=(XY3(1)-XY2(1))/L2;lanY2=(XY3(2)-XY2(2))/L2;   
G1=EA1/L1;
G2=EA2/L2;
% initial assemble global K matrix
        %  5    6
   K = [G1*lanX1^2+G2*lanX2^2];        
        f=K*u-F;
        %calculate du1 & du2
        du = Kini^-1 *(-f);
        %Calculate the residual (internal-external force deviation)
        residual = max(abs(du)); % evaluate residual
        it = it + 1; % increment iteration count
        if it == itermax
          fprintf('(-)For increment %1.0f trail iterations reached to Ultimate %1.0f\n',i,it)
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
DU1(i)=residual;I1(i)=i;IT1(i)=it;
XY3ox1(i)=[XY3(1)];XY3oy1(i)=[XY3(2)];
es1=(L1-L1i)/L1i;es2=(L2-L2i)/L2i;
ess1o1(i)=es1;ess2o1(i)=es2;
%% Force and Dispalcement for each increment
% Internal element force          
         % Displacement Transformation Matrix
        T1 = [lanX1 lanY1 0 0;
              0 0 lanX1 lanY1];
        % Stiffness Matrix for each element
        Kele1 = [G1 -G1;
                -G1 G1]; 
        Fele1 = T1'*Kele1*T1*[0;0;u(1);up];
        Fele1i1(i) = EA1*es1;%G1*[-lanX1 -lanY1 lanX1 lanY1]*[0;0;u(1);up];% Internal Force of element 1                       
        % Displacement Transformation Matrix
        T2 = [lanX2 lanY2 0 0;
              0 0 lanX2 lanY2];
        Kele2 = [G2 -G2;
                -G2 G2];        
        Fele2 = T2'*Kele2*T2*[0;0;u(1);up];
        Fele2i1(i) = EA2*es2;%G2*[-lanX2 -lanY2 lanX2 lanY2]*[0;0;u(1);up];% Internal Force of element 2                
% Force and Dispalcement for each increment
    U1N(i) = u(1);
    U2N(i) = up;
    INT_N_f1(i) = -EA1*es1*lanY1;%Fele1(1);% Base reaction of element 1
    INT_N_f2(i) = -EA2*es2*lanY2;%Fele2(1);% Base reaction of element 2
    TBSN1(i)=[INT_N_f1(i)+INT_N_f2(i)]; % Total Base Shear
    if abs(up) >= D6max;disp('      ## Displacement at [DOF (6)] reached to Ultimate Displacement ##');break;end
end
D1=[U2N'];F1=[TBSN1'];
disp('+ Large Displacement Analysis [DOF(6)] +');
disp(' X-displacement(D5)   Y-displacement(D6)');
disp('----------------------------------------')
disp([U1N' U2N'])
disp('+ ==================================== +');
disp('  Axial Load(f1)   Axial Load(f2)       ');
disp('----------------------------------------');
disp([INT_N_f1' INT_N_f2']);
disp('+======================================+');
disp('  Y-displacement(D6)   Base Reaction(kN) ');
disp('----------------------------------------');
disp([U2N' TBSN1']);
XXi1=[XY1i(1),XY3i(1),XY2i(1)];YYi1=[XY1i(2),XY3i(2),XY2i(2)];
XX1=[XY1(1),XY3(1),XY2(1)];YY1=[XY1(2),XY3(2),XY2(2)];
%% Large Displacement Analysis [DOF(5)]
disp('########################################');
disp('# Large Displacement Analysis [DOF(5)] #');
disp('########################################');
u = zeros(1,1);% initial guess value
% Gradually increase the applied displacement
for i=1:m;
up=D5*i;
XY1=[XY1i(1) XY1i(2)]; % [x y] Point 1 Coordinate
XY2=[XY2i(1) XY2i(2)]; % [x y] Point 2 Coordinate
XY3=[XY3i(1)+up XY3i(2)+u(1)]; % [x y] Point 3 Coordinate
L1 = ((XY3(1)-XY1(1))^2+(XY3(2)-XY1(2))^2)^.5; % [mm]
L2 = ((XY3i(1)-XY2i(1))^2+(XY3i(2)-XY2i(2))^2)^.5; % [mm]        
lanX1=(XY3(1)-XY1(1))/L1;lanY1=(XY3(2)-XY1(2))/L1;
lanX2=(XY3(1)-XY2(1))/L2;lanY2=(XY3(2)-XY2(2))/L2;    
G1=EA1/L1;
G2=EA2/L2;
% initial assemble global K matrix
        %  5    6
Kp = [G1*lanX1^2+G2*lanX2^2 G1*lanX1*lanY1+G2*lanX2*lanY2;
      G1*lanX1*lanY1+G2*lanX2*lanY2 G1*lanY1^2+G2*lanY2^2];
Fii=Kp(:,1)*up;
Kini = [G1*lanY1^2+G2*lanY2^2];
% Define the applied load
Fi = [P5;P6];
F=Fi-Fii;F=[F(2,1)];
        it = 0; % initialize iteration count
        residual = 100; % initialize residual
      while (residual > tolerance)
XY1=[XY1i(1) XY1i(2)]; % [x y] Point 1 Coordinate
XY2=[XY2i(1) XY2i(2)]; % [x y] Point 2 Coordinate
XY3=[XY3i(1)+up XY3i(2)+u(1)]; % [x y] Point 3 Coordinate          
L1 = ((XY3(1)-XY1(1))^2+(XY3(2)-XY1(2))^2)^.5; % [mm]
L2 = ((XY3i(1)-XY2i(1))^2+(XY3i(2)-XY2i(2))^2)^.5; % [mm]        
lanX1=(XY3(1)-XY1(1))/L1;lanY1=(XY3(2)-XY1(2))/L1;
lanX2=(XY3(1)-XY2(1))/L2;lanY2=(XY3(2)-XY2(2))/L2;   
G1=EA1/L1;
G2=EA2/L2;
% initial assemble global K matrix
        %  5    6
   K = [G1*lanY1^2+G2*lanY2^2];        
        f=K*u-F;
        %calculate du1 & du2
        du = Kini^-1 *(-f);
        %Calculate the residual (internal-external force deviation)
        residual = max(abs(du)); % evaluate residual
        it = it + 1; % increment iteration count
        if it == itermax
          fprintf('(-)For increment %1.0f trail iterations reached to Ultimate %1.0f\n',i,it)
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
DU2(i)=residual;I2(i)=i;IT2(i)=it;
XY3ox2(i)=[XY3(1)];XY3oy2(i)=[XY3(2)];
es1=(L1-L1i)/L1i;es2=(L2-L2i)/L2i;
ess1o2(i)=es1;ess2o2(i)=es2;LL1(i)=L1;LL2(i)=L2;
%% Force and Dispalcement for each increment
% Internal element force          
         % Displacement Transformation Matrix
        T1 = [lanX1 lanY1 0 0;
              0 0 lanX1 lanY1];
        % Stiffness Matrix for each element
        Kele1 = [G1 -G1;
                -G1 G1]; 
        Fele1 = T1'*Kele1*T1*[0;0;up;u(1)];
        Fele1i2(i) = EA1*es1;% Internal Force of element 1        
        % Displacement Transformation Matrix
        T2 = [lanX2 lanY2 0 0;
              0 0 lanX2 lanY2];
        Kele2 = [G2 -G2;
                -G2 G2];        
        Fele2 = T2'*Kele2*T2*[0;0;up;u(1)];
        Fele2i2(i) = EA2*es2;% Internal Force of element 2        
% Force and Dispalcement for each increment
    U1L(i) = up;
    U2L(i) = u(1);
    INT_L_f1(i) = -EA1*es1*lanX1;%Fele1(1);% Base reaction of element 1
    INT_L_f2(i) = -EA2*es2*lanX2;%Fele2(1);% Base reaction of element 2
    INT_L_f1y(i) = -EA1*es1*lanY1;INT_L_f2y(i) = -EA2*es2*lanY2;
    TBSL2(i)=[INT_L_f1(i)+INT_L_f2(i)]; % Total Base Shear
    if abs(up) >= D5max;disp('      ## Displacement at [DOF (5)] reached to Ultimate Displacement ##');break;end
end
D2=[U1L'];F2=[TBSL2'];
disp('+ Large Displacement Analysis [DOF(5)] +');
disp(' X-displacement(D4)   Y-displacement(D5)');
disp('----------------------------------------')
disp([U1L' U2L'])
disp('+======================================+');
disp('  Axial Load(f1)   Axial Load(f2)       ');
disp('----------------------------------------');
disp([INT_L_f1' INT_L_f2']);
disp('+======================================+');
disp('  Y-displacement(D5)   Base Reaction(kN) ');
disp('----------------------------------------');
disp([U1L' TBSL2']);
disp('+======================================+');
XXi20=[XY1i(1),XY3i(1),XY2i(1)];YYi20=[XY1i(2),XY3i(2),XY2i(2)];
XXi21=[XY1i(1),XY3i(1)+U1L(25),XY2i(1)];YYi21=[XY1i(2),XY3i(2),XY2i(2)];
XXi22=[XY1i(1),XY3i(1)+U1L(50),XY2i(1)];YYi22=[XY1i(2),XY3i(2),XY2i(2)];
XXi23=[XY1i(1),XY3i(1)+U1L(75),XY2i(1)];YYi23=[XY1i(2),XY3i(2),XY2i(2)];
XXi24=[XY1(1),XY3(1),XY2(1)];YYi24=[XY1(2),XY3(2),XY2(2)];
%%%  print time of computation
totaltime = cputime - starttime;
fprintf('\nTotal time (s): %7.4f \n\n',totaltime)
%% imaging
figure(1)
IMAGE=imread('PushoverLateralTruss03Ldc.jpg');
image(IMAGE);axis image;axis off;
figure(2)
p1=plot(I1,DU1,'black',I2,DU2,'--blue');grid on;set(p1,'LineWidth',2);
xlabel('increment');ylabel('Residual');
title('Residual-increment diagram','color','b');
legend('[DOF(6)]','[DOF(5)]','Location','NorthEastOutside');
figure(3)
p1=plot(I1,IT1,'black',I2,IT2,'--blue');grid on;set(p1,'LineWidth',2);
xlabel('increment');ylabel('Iteration');
title('Iteration-increment diagram','color','b');
legend('[DOF(6)]','[DOF(5)]','Location','NorthEastOutside');
figure(4)
p1=plot(lanX1N1,lanY1N1,'-.black',lanX2N1,lanY2N1,'-.blue');grid on;set(p1,'LineWidth',3);
legend('Angle of element 1','Angle of element 2','Location','NorthEastOutside');
xlabel('Angle of element during the incremental Loading [cos (a)] (rad)');ylabel('Angle of element during the incremental Loading [sin (a)] (rad)]');
title('[DOF(6)]: Angle of elements during the incremental displacement','color','b')
figure(5)
p1=plot(U1N,U2N,'oblack');grid on;set(p1,'LineWidth',3);
xlabel('Displacement-X (mm) [DOF(5)]');ylabel('Displacement-Y (mm) [DOF(6)]');
title('[DOF(6)]: Large Displacement Theory - Displacement-X and Y during the incremental displacement[DOF(6)]','color','b');
figure(6)
p1=plot(ess1o1,ess2o1,'oblack');grid on;set(p1,'LineWidth',3);
xlabel('Ele.1 Strain (mm/mm)');ylabel('Ele.2 Strain (mm/mm)');
title('Strain in element of truss during the incremental displacement [DOF(6)]  (Positive=Tension) (Negative=Compression)','color','b');
figure(7)
p1=plot(XY3ox1,XY3oy1,'oblack');grid on;set(p1,'LineWidth',3);
xlabel('Coordinate of node(3) (mm) [DOF(5)]');ylabel('Coordinate of node(3) (mm) [DOF(6)]');
title('[DOF(6)]: Large Displacement Theory - Coordinate of node(3) during the incremental displacement[DOF(6)]','color','b');
figure(8)
p1=plot(D1,F1);grid on;set(p1,'LineWidth',3);
xlabel('Displacement (mm) [DOF(6)]');ylabel('Base Reaction (kN) [DOF(2)+DOF(4)]');
title('[DOF(6)]: Base Reaction-Displacement of truss during the incremental displacement','color','b');
figure(9)
p1=plot(D1,Fele1i1,D1,Fele2i1,'r--');grid on;set(p1,'LineWidth',3);
xlabel('Displacement (mm) [DOF(6)]');ylabel('Axial Load of element (kN)');
title('Element Axial Load-Displacement of truss during the incremental displacement [DOF(6)]  (Positive=Tension) (Negative=Compression)','color','b');
legend('Analysis-ele. 1','Analysis-ele. 2','Location','NorthEastOutside');
figure(10)
p1=plot(XXi1,YYi1,XX1,YY1,'r--');grid on;set(p1,'LineWidth',3);
xlabel('X-Diemention (mm)');ylabel('Y-Diemention (mm)');
title('Shape of truss during the incremental displacement [DOF(6)]','color','b');
legend('Not loading step','Last step','Location','NorthEastOutside');
figure(11)
p1=plot(lanX1L1,lanY1L1,'-.black',lanX2L1,lanY2L1,'-.blue');grid on;set(p1,'LineWidth',3);
legend('Angle of element 1','Angle of element 2','Location','NorthEastOutside');
xlabel('Angle of element during the incremental Loading [cos (a)] (rad)');ylabel('Angle of element during the incremental Loading [sin (a)] (rad)]');
title('[DOF(5)]: Angle of elements during the incremental displacement','color','b')
figure(12)
p1=plot(U1L,U2L,'oblack');grid on;set(p1,'LineWidth',3);
xlabel('Displacement-X (mm) [DOF(5)]');ylabel('Displacement-Y (mm) [DOF(6)]');
title('[DOF(5)]: Large Displacement Theory - Displacement-X and Y during the incremental displacement','color','b');
figure(13)
p1=plot(ess1o2,ess2o2,'oblack');grid on;set(p1,'LineWidth',3);
xlabel('Ele.1 Strain (mm/mm)');ylabel('Ele.2 Strain (mm/mm)');
title('Strain in element of truss during the incremental displacement [DOF(5)]  (Positive=Tension) (Negative=Compression)','color','b');
figure(14)
p1=plot(XY3ox2,XY3oy2,'oblack');grid on;set(p1,'LineWidth',3);
xlabel('Coordinate of node(3) (mm) [DOF(5)]');ylabel('Coordinate of node(3) (mm) [DOF(6)]');
title('[DOF(5)]: Large Displacement Theory - Coordinate of node(3) during the incremental displacement[DOF(5)]','color','b');
figure(15)
p1=plot(D2,INT_L_f1y,D2,INT_L_f2y,'r--');grid on;set(p1,'LineWidth',3);
xlabel('Displacement (mm) [DOF(5)]');ylabel('Base Reaction (kN)');
title('[DOF(5)]: Base Reaction-Displacement of truss during the incremental displacement[DOF(5)]','color','b');
figure(16)
p1=plot(D2,F2);grid on;set(p1,'LineWidth',3);
xlabel('Displacement (mm) [DOF(5)]');ylabel('Base Reaction (kN) [DOF(1)+DOF(3)]');
title('[DOF(5)]: Base Reaction-Displacement of truss during the incremental displacement[DOF(5)]','color','b');
figure(17)
p1=plot(D2,Fele1i2,D2,Fele2i2,'r--');grid on;set(p1,'LineWidth',3);
xlabel('Displacement (mm) [DOF(5)]');ylabel('Axial Load of element (kN)');
title('Element Axial Load-Displacement of truss during the incremental displacement [DOF(5)]  (Positive=Tension) (Negative=Compression)','color','b');
legend('Analysis-ele. 1','Analysis-ele. 2','Location','NorthEastOutside');
figure(18)
p1=plot(XXi20,YYi20,XXi21,YYi21,'r--',XXi22,YYi22,'r--',XXi23,YYi23,'r--',XXi24,YYi24,'r--');grid on;set(p1,'LineWidth',3);
xlabel('X-Diemention (mm)');ylabel('Y-Diemention (mm)');
title('Shape of truss during the incremental displacement [DOF(5)]','color','b');
legend('Not loading step [DOF(5)=0 mm]','[DOF(5)=25 mm]','[DOF(5)=50 mm]','[DOF(5)=75 mm]','Last step [DOF(5)=100 mm]','Location','NorthEastOutside');