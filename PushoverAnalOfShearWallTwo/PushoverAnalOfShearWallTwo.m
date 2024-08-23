%***********************************************************%
%                >> IN THE NAME OF GOD <<                   %
% 1st order and 2nd order Analysis of Interaction           %
% of Shear wall with Elstic Frame with Nonlinear Semi Rigid %
% Conection Frame subjected to Pushover lateral load        %
% Small Deflection Theory                                   %
% Modified Newton Raphson Method for single increment:      %
% In this procedure Initial stiffness change in each step   %
%-----------------------------------------------------------%
%     This program is written by salar delavar ghashghaei   %
%          E-mail:salar.d.ghashghaei@gmail.com              %
%-----------------------------------------------------------%
%                                 w                         %
%                          ||||||||||||||                   %
%      Horizontal Force -> +------------+   -               %
%                          |     EIb    |   |               %
%                          |EIc         |   h               %
%                          @            @   |               %
%                          =            =   -               %
%                          |---- L -----|                   %
%***********************************************************%
close all;clear all;clc;
% Define Parameters in mm,kN
W=.005; % [kN/mm] % Distributed load Value (+ : Down)
h= 3000; % [mm] % column length
L = 6000; % [mm] % beam length
bw = 3000; % [mm] width of shear wall
Lw = 3000; % [mm] Length of shear wall
Ew = 200; % [kN/mm^2] modulus of elasticity of shear wall
EIc = 200*100^4/12; % [kN.mm^2] column
EAc = 200*10000; % [kN]
EIb = 200*50^4/12; % [kN.mm^2] beam
EAb = 200*(50)^2; % [kN]
m = 2000; % number of calculation (Load Steps)
itermax = 400;% maximum number of iterations
tolerance = 1e-9; % specified tolerance for convergence
u = zeros(9,1);% initial guess value
P3 =0; % [kN.mm] Moment [DOF (3)]
P6 =0; % [kN.mm] Moment [DOF (6)]
P7=50; % [kN] Horizental Force [DOF (7)] Incremantal Loading
P8=-.5*W*L; % [kN] Vertical Force [DOF (8)]
P9 =-(W*L^2)/12; % [kN.mm] Moment [DOF (9)]
P10=0; % [kN] Horizental Force [DOF (10)]
P11=-.5*W*L; % [kN] Vertical Force [DOF (11)]
P12=+(W*L^2)/12; % [kN.mm] Moment [DOF (12)]
lanXc=0;lanYc=1;
lanXb=1;lanYb=0;
%% Nonlinear Rotational Spring of columns
tyc=.08; % Yield rotaion
Myc=40e+3; % Yield moment
tuc=.25;  % Ultimate rotation
Muc=1.25*Myc;  % Ultimate moment
nc = 5; % Moment-rotation shape parameter
Rkic=Myc/tyc;
Rkpc=(Muc-Myc)/(tuc-tyc);
%% Nonlinear Displacement Spring of Shear Wall - Flexural Behavior
hw=h;Iw= bw*(Lw^3)/12;As=(2*Iw)/Lw^2;aw=(2*Iw)/(bw*Lw^3);
dywf=15; % Yield displacement
Fywf=20e+3; % Yield Shear Force
duwf=30;  % Ultimate displacement
Fuwf=25e+3;  % Ultimate Shear Force
nwf = 5; % Shear Force-displacement shape parameter
Rkiwf=Fywf/dywf;
Rkpwf=(Fuwf-Fywf)/(duwf-dywf);
%% Nonlinear Displacement Spring of Shear Wall - Shear Behavior
hc = hw/3; % [mm] 
dyws=5; % Yield displacement
Fyws=20e+3; % Yield Shear Force
duws=50;  % Ultimate displacement
Fuws=30e+3;  % Ultimate Shear Force
nws = 5; % Shear Force-displacement shape parameter
Rkiws=Fyws/dyws;
Rkpws=(Fuws-Fyws)/(duws-dyws);
% Notice: Semi-Rigid Spring for column Connection Number 1 
%% monitor cpu time
starttime = cputime;
%% First-order Nonlinear Analysis
disp('#################################################');
disp('#         First-order Nonlinear Analysis        #');
disp('#################################################');
% Element stiffness coefficient for element 1 and 2
A=4*EIc/h;B=6*EIc/h^2;C=2*EIc/h;D=12*EIc/h^3;G=EAc/h;
% Element stiffness coefficient for element 3
AA=4*EIb/L;BB=6*EIb/L^2;CC=2*EIb/L;DD=12*EIb/L^3;GG=EAb/L;
% Shear wall stiffness coefficient for element 3
AA=4*EIb/L;BB=6*EIb/L^2;CC=2*EIb/L;DD=12*EIb/L^3;GG=EAb/L;
for i=1:m
a1=((Rkic-Rkpc)/((1+(abs((Rkic*u(1))/Myc))^nc)^(1/nc)))+Rkpc;% Column hinge [DOF (3)]
a2=((Rkic-Rkpc)/((1+(abs((Rkic*u(2))/Myc))^nc)^(1/nc)))+Rkpc;% Column hinge [DOF (6)]
a5=((Rkiwf-Rkpwf)/((1+(abs((Rkiwf*u(9))/Fywf))^nwf)^(1/nwf)))+Rkpwf;% Shear-wall Stiffness - Flexural Behavior [DOF (10)]
a6=((Rkiws-Rkpws)/((1+(abs((Rkiws*u(6))/Fyws))^nws)^(1/nws)))+Rkpws;% Shear-wall Stiffness - Shear Behavior [DOF (13)]
% initial assemble global K matrix
        %  3    6   7    8    9    10     11    12   13
  Kinit = [A+a1 0 B*lanYc -B*lanXc C 0 0 0 0;
           0 A+a2 0 0 0 B*lanYc -B*lanXc C 0;
           B*lanYc 0 (G*lanXc^2+D*lanYc^2)+(GG*lanXb^2+DD*lanYb^2) (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb B*lanYc -(GG*lanXb^2+DD*lanYb^2) -(GG-DD)*lanXb*lanYb 0 0;
          -B*lanXc 0 (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb (G*lanYc^2+D*lanXc^2)+(GG*lanYb^2+DD*lanXb^2) -B*lanXc -(GG-DD)*lanXb*lanYb -(GG*lanYb^2+DD*lanXb^2) 0 0;
           C 0 B*lanYc -B*lanXc A+AA BB*lanYb -BB*lanXb CC 0;
           0 B*lanYc -(GG*lanXb^2+DD*lanYb^2) -(GG-DD)*lanXb*lanYb BB*lanYb (G*lanXc^2+D*lanYc^2)+(GG*lanXb^2+DD*lanYb^2)+a6 (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb B*lanYc -a6*(hw-hc);
           0 -B*lanXc -(GG-DD)*lanXb*lanYb -(GG*lanYb^2+DD*lanXb^2) -BB*lanXb (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb (G*lanYc^2+D*lanXc^2)+(GG*lanYb^2+DD*lanXb^2) -B*lanXc 0;
           0 C 0 0 0 B*lanYc -B*lanXc A+AA 0;
           0 0 0 0 CC -a6*(hw-hc) 0 0 a6*(hw-hc)^2+a5*.5*Lw^2];
% Define the applied load
F = [P3;P6;P7*i;P8;P9;P10;P11;P12;0];
it = 0; % initialize iteration count
residual = 100; % initialize residual
% calculate Force
while (residual > tolerance)
% assemble global K matrix
a1=((Rkic-Rkpc)/((1+(abs((Rkic*u(1))/Myc))^nc)^(1/nc)))+Rkpc;
a2=((Rkic-Rkpc)/((1+(abs((Rkic*u(2))/Myc))^nc)^(1/nc)))+Rkpc;
a5=((Rkiwf-Rkpwf)/((1+(abs((Rkiwf*u(9))/Fywf))^nwf)^(1/nwf)))+Rkpwf;
a6=((Rkiws-Rkpws)/((1+(abs((Rkiws*u(6))/Fyws))^nws)^(1/nws)))+Rkpws;
      K = [A+a1 0 B*lanYc -B*lanXc C 0 0 0 0;
           0 A+a2 0 0 0 B*lanYc -B*lanXc C 0;
           B*lanYc 0 (G*lanXc^2+D*lanYc^2)+(GG*lanXb^2+DD*lanYb^2) (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb B*lanYc -(GG*lanXb^2+DD*lanYb^2) -(GG-DD)*lanXb*lanYb 0 0;
          -B*lanXc 0 (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb (G*lanYc^2+D*lanXc^2)+(GG*lanYb^2+DD*lanXb^2) -B*lanXc -(GG-DD)*lanXb*lanYb -(GG*lanYb^2+DD*lanXb^2) 0 0;
           C 0 B*lanYc -B*lanXc A+AA BB*lanYb -BB*lanXb CC 0;
           0 B*lanYc -(GG*lanXb^2+DD*lanYb^2) -(GG-DD)*lanXb*lanYb BB*lanYb (G*lanXc^2+D*lanYc^2)+(GG*lanXb^2+DD*lanYb^2)+a6 (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb B*lanYc -a6*(hw-hc);
           0 -B*lanXc -(GG-DD)*lanXb*lanYb -(GG*lanYb^2+DD*lanXb^2) -BB*lanXb (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb (G*lanYc^2+D*lanXc^2)+(GG*lanYb^2+DD*lanXb^2) -B*lanXc 0;
           0 C 0 0 0 B*lanYc -B*lanXc A+AA 0;
           0 0 0 0 CC -a6*(hw-hc) 0 0 a6*(hw-hc)^2+a5*.5*Lw^2];
        f=F-K*u;
        %calculate du
        du = Kinit^-1 *(f);
        u = u+du; % update u
        %Calculate the residual (internal-external force deviation)
        residual = max(abs(du)); % evaluate residual
        it = it + 1; % increment iteration count
        if it == itermax
          fprintf('(-)For increment %1.0f trail iteration reached to Ultimate %1.0f\n',i,it)
             disp('    ## The solution for this step is not converged ##')
        end
end      
%% Force and Dispalcement for each increment
    F7(i) = F(3);
    U3(i) = u(1);
    U6(i) = u(2);
    U7(i) = u(3);
    U8(i) = u(4);
    U9(i) = u(5);
    U10(i) = u(6);
    U11(i) = u(7);
    U12(i) = u(8);
    DU1(i)=residual;I1(i)=i;IT1(i)=it;
    % Moment-Rotation of spring during the analysis
    uc1(i)=u(1);Hc1(i)=(((Rkic-Rkpc)*uc1(i))/((1+(abs((Rkic*uc1(i))/Myc))^nc)^(1/nc)))+Rkpc*uc1(i);
    uc2(i)=u(2);Hc2(i)=(((Rkic-Rkpc)*uc2(i))/((1+(abs((Rkic*uc2(i))/Myc))^nc)^(1/nc)))+Rkpc*uc2(i);
    uw13(i)=u(9);Hw13(i)=(((Rkiwf-Rkpwf)*uw13(i))/((1+(abs((Rkiwf*uw13(i))/Fywf))^nwf)^(1/nwf)))+Rkpwf*uw13(i);
    uw7(i)=u(6);Hw7(i)=(((Rkiws-Rkpws)*uw7(i))/((1+(abs((Rkiws*uw7(i))/Fyws))^nws)^(1/nws)))+Rkpws*uw7(i);
   if it < itermax
   fprintf('(+)It is converged in %1.0f iterations for increment %1.0f\n',it,i)
   end
      if 1 == isnan(u) ;break;end % Check for Softening
        if abs(u(1)) >= tuc
        disp('      ## spring at support reached to Ultimate Rotation ##')
        break
        end
end
D1=[0;U7'];F1=[0;F7'];
% Structural Lateral Stiffness [DOF(7)]
s=size(D1,1);for i=1:s-1;Kz1(i)=abs((F1(i))/(D1(i)));end
% Coordinate of element
corx=[0 0 L L];cory=[0 h h 0];
% displaying shape of deflection along the length for First-order Nonlinear Analysis
% For element 1
for j=1:11;
    y1(j)=((h+u(4))/10)*(j-1);
    N1 =(1/h^3)*(h^3-3*h*y1(j)^2+2*y1(j)^3);
    N2 =(1/h^2)*(y1(j)*h^2-2*h*y1(j)^2+y1(j)^3);
    N3 =(1/h^3)*(3*h*y1(j)^2-2*y1(j)^3);
    N4 =(1/h^2)*(-h*y1(j)^2+y1(j)^3);
    x1(j) =[N1 N2 N3 N4]*[0;u(1);u(3);u(5)];
end
% For element 3
for j=1:11;
    x2(j)=(L/10)*(j-1)+u(3);
    N1 =(1/L^3)*(L^3-3*L*x2(j)^2+2*x2(j)^3);
    N2 =(1/L^2)*(x2(j)*L^2-2*L*x2(j)^2+x2(j)^3);
    N3 =(1/L^3)*(3*L*x2(j)^2-2*x2(j)^3);
    N4 =(1/L^2)*(-L*x2(j)^2+x2(j)^3);
    y2(j) =[N1 N2 N3 N4]*[u(4);u(5);u(7);u(8)]+h;
end
% For element 2
for j=1:11;
    y3(j)=((h+u(7))/10)*(j-1);
    N1 =(1/h^3)*(h^3-3*h*y3(j)^2+2*y3(j)^3);
    N2 =(1/h^2)*(y3(j)*h^2-2*h*y3(j)^2+y3(j)^3);
    N3 =(1/h^3)*(3*h*y3(j)^2-2*y3(j)^3);
    N4 =(1/h^2)*(-h*y3(j)^2+y3(j)^3);
    x3(j) =[N1 N2 N3 N4]*[0;u(2);u(6);u(8)]+L;
end
X1N1=x1;Y1N1=y1;X2N1=x2;Y2N1=y2;X3N1=x3;Y3N1=y3;
%% Second order Nonlinear Analysis
disp('#################################################');
disp('#        Second-order Nonlinear Analysis        #');
disp('#################################################');
U = zeros(9,1);% initial guess value
% Element stiffness coefficient for element 1 and 2
A=(4*EIc/h)+(2*P8*h/15);B=(6*EIc/h^2)+(P8/10);C=(2*EIc/h)-(P8*h/30);D=(12*EIc/h^3)+(6*P8/(5*h));G=(EAc/h)+(P8/h);
% Element stiffness coefficient for element 3
AA=(4*EIb/L);BB=(6*EIb/L^2);CC=(2*EIb/L);DD=(12*EIb/h^3);GG=(EAb/L);
for i=1:m
a1=((Rkic-Rkpc)/((1+(abs((Rkic*U(1))/Myc))^nc)^(1/nc)))+Rkpc;
a2=((Rkic-Rkpc)/((1+(abs((Rkic*U(2))/Myc))^nc)^(1/nc)))+Rkpc;
a5=((Rkiwf-Rkpwf)/((1+(abs((Rkiwf*U(9))/Fywf))^nwf)^(1/nwf)))+Rkpwf;
a6=((Rkiws-Rkpws)/((1+(abs((Rkiws*U(6))/Fyws))^nws)^(1/nws)))+Rkpws;
% initial assemble global K matrix
        %  3    6   7    8    9    10     11    12   13
  Kinit = [A+a1 0 B*lanYc -B*lanXc C 0 0 0 0;
           0 A+a2 0 0 0 B*lanYc -B*lanXc C 0;
           B*lanYc 0 (G*lanXc^2+D*lanYc^2)+(GG*lanXb^2+DD*lanYb^2) (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb B*lanYc -(GG*lanXb^2+DD*lanYb^2) -(GG-DD)*lanXb*lanYb 0 0;
          -B*lanXc 0 (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb (G*lanYc^2+D*lanXc^2)+(GG*lanYb^2+DD*lanXb^2) -B*lanXc -(GG-DD)*lanXb*lanYb -(GG*lanYb^2+DD*lanXb^2) 0 0;
           C 0 B*lanYc -B*lanXc A+AA BB*lanYb -BB*lanXb CC 0;
           0 B*lanYc -(GG*lanXb^2+DD*lanYb^2) -(GG-DD)*lanXb*lanYb BB*lanYb (G*lanXc^2+D*lanYc^2)+(GG*lanXb^2+DD*lanYb^2)+a6 (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb B*lanYc -a6*(hw-hc);
           0 -B*lanXc -(GG-DD)*lanXb*lanYb -(GG*lanYb^2+DD*lanXb^2) -BB*lanXb (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb (G*lanYc^2+D*lanXc^2)+(GG*lanYb^2+DD*lanXb^2) -B*lanXc 0;
           0 C 0 0 0 B*lanYc -B*lanXc A+AA 0;
           0 0 0 0 CC -a6*(hw-hc) 0 0 a6*(hw-hc)^2+a5*.5*Lw^2];
% Define the applied load
F = [P3;P6;P7*i;P8;P9;P10;P11;P12;0];
it = 0; % initialize iteration count
residual = 100; % initialize residual
% calculate Force
while (residual > tolerance)
% assemble global K matrix
a1=((Rkic-Rkpc)/((1+(abs((Rkic*U(1))/Myc))^nc)^(1/nc)))+Rkpc;
a2=((Rkic-Rkpc)/((1+(abs((Rkic*U(2))/Myc))^nc)^(1/nc)))+Rkpc;
a5=((Rkiwf-Rkpwf)/((1+(abs((Rkiwf*U(9))/Fywf))^nwf)^(1/nwf)))+Rkpwf;
a6=((Rkiws-Rkpws)/((1+(abs((Rkiws*U(6))/Fyws))^nws)^(1/nws)))+Rkpws;
      K = [A+a1 0 B*lanYc -B*lanXc C 0 0 0 0;
           0 A+a2 0 0 0 B*lanYc -B*lanXc C 0;
           B*lanYc 0 (G*lanXc^2+D*lanYc^2)+(GG*lanXb^2+DD*lanYb^2) (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb B*lanYc -(GG*lanXb^2+DD*lanYb^2) -(GG-DD)*lanXb*lanYb 0 0;
          -B*lanXc 0 (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb (G*lanYc^2+D*lanXc^2)+(GG*lanYb^2+DD*lanXb^2) -B*lanXc -(GG-DD)*lanXb*lanYb -(GG*lanYb^2+DD*lanXb^2) 0 0;
           C 0 B*lanYc -B*lanXc A+AA BB*lanYb -BB*lanXb CC 0;
           0 B*lanYc -(GG*lanXb^2+DD*lanYb^2) -(GG-DD)*lanXb*lanYb BB*lanYb (G*lanXc^2+D*lanYc^2)+(GG*lanXb^2+DD*lanYb^2)+a6 (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb B*lanYc -a6*(hw-hc);
           0 -B*lanXc -(GG-DD)*lanXb*lanYb -(GG*lanYb^2+DD*lanXb^2) -BB*lanXb (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb (G*lanYc^2+D*lanXc^2)+(GG*lanYb^2+DD*lanXb^2) -B*lanXc 0;
           0 C 0 0 0 B*lanYc -B*lanXc A+AA 0;
           0 0 0 0 CC -a6*(hw-hc) 0 0 a6*(hw-hc)^2+a5*.5*Lw^2];
        f=F-K*U;
        %calculate du
        du = Kinit^-1 *(f);
        U = U+du; % update u
        %Calculate the residual (internal-external force deviation)
        residual = max(abs(du)); % evaluate residual
        it = it + 1; % increment iteration count
        if it == itermax
          fprintf('(-)For increment %1.0f trail iteration reached to Ultimate %1.0f\n',i,it)
             disp('    ## The solution for this step is not converged ##')
        end
end
%% Force and Dispalcement for each increment
    FI7(i) = F(3);
    UI3(i) = U(1);
    UI6(i) = U(2);
    UI7(i) = U(3);
    UI8(i) = U(4);
    UI9(i) = U(5);
    UI10(i) = U(6);
    UI11(i) = U(7);
    UI12(i) = U(8);
    DU2(i)=residual;I2(i)=i;IT2(i)=it;
    % Moment-Rotation of spring during the analysis
    Uc1(i)=U(1);HC1(i)=(((Rkic-Rkpc)*Uc1(i))/((1+(abs((Rkic*Uc1(i))/Myc))^nc)^(1/nc)))+Rkpc*Uc1(i);
    Uc2(i)=U(2);HC2(i)=(((Rkic-Rkpc)*Uc2(i))/((1+(abs((Rkic*Uc2(i))/Myc))^nc)^(1/nc)))+Rkpc*Uc2(i);
    Uw13(i)=U(9);HW13(i)=(((Rkiwf-Rkpwf)*Uw13(i))/((1+(abs((Rkiwf*Uw13(i))/Fywf))^nwf)^(1/nwf)))+Rkpwf*Uw13(i);
    Uw7(i)=U(6);HW7(i)=(((Rkiws-Rkpws)*Uw7(i))/((1+(abs((Rkiws*Uw7(i))/Fyws))^nws)^(1/nws)))+Rkpws*Uw7(i);
   if it < itermax
   fprintf('(+)It is converged in %1.0f iterations for increment %1.0f\n',it,i)
   end
   if 1 == isnan(U) ;break;end % Check for Softening
        if abs(U(1)) >= tuc
        disp('      ## spring at support reached to Ultimate Rotation ##')
        break
        end
end
D2=[0;UI7'];F2=[0;FI7'];
% Structural Lateral Stiffness [DOF(7)]
s=size(D2,1);for i=1:s-1;Kz2(i)=abs((F2(i))/(D2(i)));end
% displaying shape of deflection along the length for First-order Nonlinear Analysis
% For element 1
for j=1:11;
    y1(j)=((h+U(4))/10)*(j-1);
    N1 =(1/h^3)*(h^3-3*h*y1(j)^2+2*y1(j)^3);
    N2 =(1/h^2)*(y1(j)*h^2-2*h*y1(j)^2+y1(j)^3);
    N3 =(1/h^3)*(3*h*y1(j)^2-2*y1(j)^3);
    N4 =(1/h^2)*(-h*y1(j)^2+y1(j)^3);
    x1(j) =[N1 N2 N3 N4]*[0;U(1);U(3);U(5)];
end
% For element 3
for j=1:11;
    x2(j)=(L/10)*(j-1)+U(3);
    N1 =(1/L^3)*(L^3-3*L*x2(j)^2+2*x2(j)^3);
    N2 =(1/L^2)*(x2(j)*L^2-2*L*x2(j)^2+x2(j)^3);
    N3 =(1/L^3)*(3*L*x2(j)^2-2*x2(j)^3);
    N4 =(1/L^2)*(-L*x2(j)^2+x2(j)^3);
    y2(j) =[N1 N2 N3 N4]*[U(4);U(5);U(7);U(8)]+h;
end
% For element 2
for j=1:11;
    y3(j)=((h+U(7))/10)*(j-1);
    N1 =(1/h^3)*(h^3-3*h*y3(j)^2+2*y3(j)^3);
    N2 =(1/h^2)*(y3(j)*h^2-2*h*y3(j)^2+y3(j)^3);
    N3 =(1/h^3)*(3*h*y3(j)^2-2*y3(j)^3);
    N4 =(1/h^2)*(-h*y3(j)^2+y3(j)^3);
    x3(j) =[N1 N2 N3 N4]*[0;U(2);U(6);U(8)]+L;
end
X1N2=x1;Y1N2=y1;X2N2=x2;Y2N2=y2;X3N2=x3;Y3N2=y3;
%% First-order Nonlinear Analysis bilinear fitting
SIZE=size(D1,1);
for i=1:SIZE-1;
    hh(i) = D1(i+1)-D1(i);
    Aa(i)=(F1(i)+F1(i+1))*0.5*hh(i);
end
Area=sum(Aa);k0 =F1(500)/D1(500);
fiy = (F1(i+1)*max(D1)*0.5-Area)/(F1(i+1)*0.5 - k0*max(D1)*0.5);
Fy = k0*fiy;
X1 = [0 fiy max(D1)];Y1 = [0 Fy F1(i+1)];
%% Second-order Nonlinear Analysis bilinear fitting
SIZE=size(D2,1);
for i=1:SIZE-1;
    hhl(i) = D2(i+1)-D2(i);
    Aal(i)=(F2(i)+F2(i+1))*0.5*hhl(i);
end
Area=sum(Aal);k0 =F2(500)/D2(500);
fiy = (F2(i+1)*max(D2)*0.5-Area)/(F2(i+1)*0.5 - k0*max(D2)*0.5);
Fy = k0*fiy;
X2 = [0 fiy max(D2)];Y2 = [0 Fy F2(i+1)];
%% Column - Semi Rigid Hinge Ploting
for i=1:20;
    Uc=tuc/20;uc(i)=Uc*i;
    Hc(i)=(((Rkic-Rkpc)*uc(i))/((1+(abs((Rkic*uc(i))/Myc))^nc)^(1/nc)))+Rkpc*uc(i);
end
uc =[0 uc];Hc =[0 Hc];
% Semi Rigid Column Hinge Section Ductility Rito
SIZE=size(uc,2);
for i=1:SIZE-1;
    hhg(i) = uc(i+1)-uc(i);
    Aag(i)=(Hc(i)+Hc(i+1))*0.5*hhg(i);
end
Area=sum(Aag);k0 =Hc(5)/uc(5);
fiy = (Hc(i+1)*max(uc)*0.5-Area)/(Hc(i+1)*0.5 - k0*max(uc)*0.5);
Fy = k0*fiy;
Xuc = [0 fiy max(uc)];YHc = [0 Fy Hc(i+1)];
%% Shear Wall - Hinge Ploting - Displacement
for i=1:20;
    Uc=duws/20;uws(i)=Uc*i;
    Hws(i)=(((Rkiws-Rkpws)*uws(i))/((1+(abs((Rkiws*uws(i))/Fyws))^nws)^(1/nws)))+Rkpws*uws(i);
end
uws =[0 uws];Hws =[0 Hws];
% Shear Wall - Hinge Ploting - Rotation
for i=1:20;
    Uc=duwf/20;uwf(i)=Uc*i;
    Hwf(i)=(((Rkiwf-Rkpwf)*uwf(i))/((1+(abs((Rkiwf*uwf(i))/Fywf))^nwf)^(1/nwf)))+Rkpwf*uwf(i);
end
uwf =[0 uwf];Hwf =[0 Hwf];
%% Ductility Rito if Semi-Rigid Column Connection
Z1=Xuc(3)/Xuc(2);
%% Ductility Rito if Frame
Z3=X1(3)/X1(2);Z4=X2(3)/X2(2);
%% Over Strength Ratio of Frame
Z5=Y1(3)/Y1(2);Z6=Y2(3)/Y2(2);
%% Structural Stiffness changing During the analysis
Kyf=Y1(2)/X1(2);Ktf=(Y1(3)-Y1(2))/(X1(3)-X1(2));
Kys=Y2(2)/X2(2);Kts=(Y2(3)-Y2(2))/(X2(3)-X2(2));
fprintf('======================= Result ========================\n');
fprintf('=== 1st-order Nonlinear ==+== 2nd-order Nonlinear =====\n');
fprintf('Disp.(D7) Base Shear(D1+D4) Disp.(D7) Base Shear(D1+D4)\n');
fprintf('=======================================================\n');
fprintf('  (mm)         (kN)          (mm)          (kN)\n');
fprintf('-------------------------------------------------------\n');
disp([X1' Y1' X2' Y2'])
fprintf('=======================================================\n');
fprintf('\n')
disp('+----------------------------------------------------------------------+')
fprintf(' Semi-Rigid Column Connection Ductility Rito is: %6.3f\n',Z1)
fprintf(' 1st-order Nonlinear Ductility Rito is (Du/Dy): %6.3f\n',Z3)
fprintf(' 2nd-order Nonlinear Ductility Rito is (Du/Dy): %5.3f\n',Z4)
fprintf(' 1st-order Nonlinear Over Strength Ratio is (Fu/Fy): %6.3f\n',Z5)
fprintf(' 2nd-order Nonlinear Over Strength Ratio is (Fu/Fy): %5.3f\n',Z6)
fprintf(' 1st-order Nonlinear Initial Strucural stiffness is (Ke): %6.3f [kN/mm]\n',Kyf)
fprintf(' 1st-order Nonlinear Tangent Strucural stiffness is (Kt): %5.3f [kN/mm]\n',Ktf)
fprintf(' 2nd-order Nonlinear Initial Strucural stiffness is (Ke): %6.3f [kN/mm]\n',Kys)
fprintf(' 2nd-order Nonlinear Tangent Strucural stiffness is (Kt): %5.3f [kN/mm]\n',Kts)
disp('+----------------------------------------------------------------------+')
%% Plot
figure(1)
IMAGE=imread('PushoverAnalOfNonlinSemiRigidSpringWithElaFrameUnderDistriLoadSTIFFNESS1.jpg');
image(IMAGE);axis image;axis off;
 title('Stiffness Matrix for element 1 & 2 - Columns','FontSize',12,'color','b')
figure(2)
IMAGE=imread('PushoverAnalTwoSTIFFNESS2.jpg');
image(IMAGE);axis image;axis off;
 title('Stiffness Matrix for element 3 - Beam','FontSize',12,'color','b')
figure(3)
IMAGE=imread('PushoverAnalOfShearWallTwo.jpg');
image(IMAGE);axis image;axis off;
figure(4)
IMAGE=imread('PushoverAnalOfShearWallNonlinSemiRigSpriUnderDistriLoadl.jpg');
image(IMAGE);axis image;axis off;
figure(5)
IMAGE=imread('Modified Newton Raphson Method for single increment.jpg');
image(IMAGE);axis image;axis off;
figure(6)
p1=plot(I1,DU1,'black',I2,DU2,'--blue');grid on;set(p1,'LineWidth',2);
xlabel('increment');ylabel('Residual');
title('Residual-increment diagram','color','b');
legend('1st-order Nonlinear Analysis','2nd-order Nonlinear Analysis','Location','NorthEastOutside');
figure(7)
p1=plot(I1,IT1,'black',I2,IT2,'--blue');grid on;set(p1,'LineWidth',2);
xlabel('increment');ylabel('Iteration');
title('Iteration-increment diagram','color','b');
legend('1st-order Nonlinear Analysis','2nd-order Nonlinear Analysis','Location','NorthEastOutside');
figure(8)
p1=plot(uc,Hc,'b',Xuc,YHc,'--r');grid on;set(p1,'LineWidth',3);
legend('Column Hinge','Column Hinge Curve bilinear fitted','Location','NorthEastOutside');
xlabel('Rotation (rad)');ylabel('Moment (kN.mm)');
title('Moment-Rotation Relation of Nonlin. Semi-Rigid Connection','color','b')
figure(9)
p1=plot(uws,Hws,'--b',uwf,Hwf,'--g');grid on;set(p1,'LineWidth',3);
legend('Shear wall-Shear behavior','Shear wall-flextural behavior','Location','NorthEastOutside');
xlabel('Displacement (mm)');ylabel('Force (kN)');
title('Nonlin. behavior of Force-Displacement Relation','color','b')
figure(10)
p1=plot([0 uc1],[0 Hc1],[0 uc2],[0 Hc2],[0 Uc1],[0 HC1],'--',[0 Uc2],[0 HC2],'--');grid on;set(p1,'LineWidth',3);
legend('1st-order anal. Hing DOF(3)','1st-order anal. Hing DOF(6)','2nd-order anal. Hing DOF(3)','2nd-order anal. Hing DOF(6)','Location','NorthEastOutside');
xlabel('Rotation (rad)');ylabel('Moment (kN.mm)');
title('Moment-Rotation diagram of Connection During the Analysis','color','b');
figure(11)
p1=plot([0 uw13],[0 Hw13],'b',[0 uw7],[0 Hw7],'y',[0 Uw13],[0 HW13],'--g',[0 Uw7],[0 HW7],'--r');grid on;set(p1,'LineWidth',3);
legend('1st-order Shear wall-Shear behavior','1st-order Shear wall-flextural behavior','2nd-order Shear wall-Shear behavior','2nd-order Shear wall-flextural behavior','Location','NorthEastOutside');
xlabel('Displacement (mm)');ylabel('Force (kN)');
title('Force-Displacement Relation of Nonlin. behavior of Shear-Wall During the analysis','color','b')
figure(12)
p1=plot(D1,F1,'b',X1,Y1,'r--',D2,F2,'g',X2,Y2,'--');grid on;set(p1,'LineWidth',3);
legend('1st-order Nonlinear Analysis','1st-order Nonlinear Curve fitted','2nd-order Nonlinear Analysis','2nd-order Nonlinear Curve fitted','Location','NorthEastOutside');
xlabel('Displacement (mm) [DOF (7)]');ylabel('Force (kN) [DOF (7)]');
title(['Force-Disp. diagram of Frame with Nonlin. Semi-Rigid Connection',' 1-order Nonlin. Ductility Rito: ',num2str(Z3),'   2-order Nonlin. Ductility Rito: ',num2str(Z4)],'color','b')
figure(13)
p1=plot(U7,Kz1,'black',UI7,Kz2,'b--');grid on;set(p1,'LineWidth',3);
legend('1st-order Nonlinear Analysis','2nd-order Nonlinear Analysis','Location','NorthEastOutside');
xlabel('Displacement (mm) [DOF (7)]');ylabel('Structural Lateral Stiffness (kN/mm) [DOF (7)]');
title('Structural Lateral Stiffness-Displacement diagram of Elastic Frame with Nonlin. Semi-Rigid Connection','color','b')
figure(14)
p1=plot(corx,cory,'--b',X1N1,Y1N1,'--r',X2N1,Y2N1,'--r',X3N1,Y3N1,'--r',X1N2,Y1N2,'-.g',X2N2,Y2N2,'-.g',X3N2,Y3N2,'-.g');grid on;set(p1,'LineWidth',3);
hold on
fill([L,L+U(6),L+Lw+U(6),L+Lw], [0,hw,hw,0], 'b','Parent', gca);
text(L+Lw*.5,hw*.5,'Shear Wall','color','w');
hold off
legend('BLUE: Not Loading','RED: 1st-order Nonlinear Analysis','GREEN: 2nd-order Nonlinear Analysis','BLUE BOX: Shear-Wall','Location','NorthOutside');
xlabel('X (mm)');ylabel('Y (mm)');
title('Small Deflection Theory - Last Step Displacement of Nonlin. behavior of Shear-Wall and Elastic Frame with Nonlin. Semi-Rigid Connection','color','b')
%%  print time of computation
totaltime = cputime - starttime;
fprintf('\nTotal time (s)= %7.4f \n\n',totaltime)
