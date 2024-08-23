%***********************************************************%
%                >> IN THE NAME OF GOD <<                   %
% 1st order and 2nd order Analysis of Braced with Elastic   %
% Frame with Nonlinear Semi Rigid Connection Frame          %
%  subjected to Pushover lateral load                       %
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
%                          |  \      /  |   |               %
%                          |      X     |   h               %
%                          @  /      \  @   |               %
%                          =            =   -               %
%                          |---- L -----|                   %
%***********************************************************%
close all;clear all;clc;
% Define Parameters in mm,kN
W=.03; % [kN/mm] % Distributed load Value (+ : Down)
h= 3000; % [mm] % column length
L = 6000; % [mm] % beam length
EIc = 200*100^4/12; % [kN.mm^2] column
EAc = 200*10000; % [kN]
EIb = 200*50^4/12; % [kN.mm^2] beam
EAb = 200*(50)^2; % [kN]
m = 3000; % number of calculation (Load Steps)
itermax = 400;% maximum number of iterations
tolerance = 1e-9; % specified tolerance for convergence
u = zeros(8,1);% initial guess value
P3 =0; % [kN.mm] Moment [DOF (3)]
P6 =0; % [kN.mm] Moment [DOF (6)]
P7=2; % [kN] Horizental Force [DOF (7)] Incremantal Loading
P8=-.5*W*L; % [kN] Vertical Force [DOF (8)]
P9 =-(W*L^2)/12; % [kN.mm] Moment [DOF (9)]
P10=0; % [kN] Horizental Force [DOF (10)]
P11=-.5*W*L; % [kN] Vertical Force [DOF (11)]
P12=+(W*L^2)/12; % [kN.mm] Moment [DOF (12)]
lanXc=0;lanYc=1;% Column
lanXb=1;lanYb=0;% beam
lanX4=L/sqrt(L^2+h^2);lanY4=h/sqrt(L^2+h^2);% brace (4)
lanX5=-L/sqrt(L^2+h^2);lanY5=h/sqrt(L^2+h^2);% brace (4)
%% Nonlinear Rotational Spring of columns
tyc=.08; % Yield rotaion
Myc=40e+3; % Yield moment
tuc=.25;  % Ultimate rotation
Muc=1.25*Myc;  % Ultimate moment
nc = 9; % Moment-rotation shape parameter
Rkic=Myc/tyc;
Rkpc=(Muc-Myc)/(tuc-tyc);
%% Nonlinear Axial displacement Spring of Brace
dyb=10; % Yield displacement [mm]
Fyb=1000; % Yield Shear Force [kN]
dub=35;  % Ultimate displacement [mm]
Fub=1.25*Fyb;  % Ultimate Shear Force [kN]
n = 9; % Shear Force-displacement shape parameter
Rki=Fyb/dyb;
Rkp=(Fub-Fyb)/(dub-dyb);
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
for i=1:m
a1=((Rkic-Rkpc)/((1+(abs((Rkic*u(1))/Myc))^nc)^(1/nc)))+Rkpc;% Column hinge [DOF (3)]
a2=((Rkic-Rkpc)/((1+(abs((Rkic*u(2))/Myc))^nc)^(1/nc)))+Rkpc;% Column hinge [DOF (6)]
a5=((Rki-Rkp)/((1+(abs((Rki*u(6))/Fyb))^n)^(1/n)))+Rkp;% Brace stiffness [DOF (10)]
a6=((Rki-Rkp)/((1+(abs((Rki*u(7))/Fyb))^n)^(1/n)))+Rkp;% Brace stiffness [DOF (11)]
a7=((Rki-Rkp)/((1+(abs((Rki*u(3))/Fyb))^n)^(1/n)))+Rkp;% Brace stiffness [DOF (7)]
a8=((Rki-Rkp)/((1+(abs((Rki*u(4))/Fyb))^n)^(1/n)))+Rkp;% Brace stiffness [DOF (8)]
% initial assemble global K matrix
        %  3    6   7    8    9    10     11    12
  Kinit = [A+a1 0 B*lanYc -B*lanXc C 0 0 0;
           0 A+a2 0 0 0 B*lanYc -B*lanXc C;
           B*lanYc 0 (G*lanXc^2+D*lanYc^2)+(GG*lanXb^2+DD*lanYb^2)+a7*lanX5^2 (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb B*lanYc+a8*lanX5*lanY5 -(GG*lanXb^2+DD*lanYb^2) -(GG-DD)*lanXb*lanYb 0;
          -B*lanXc 0 (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb+a7*lanX5*lanY5 (G*lanYc^2+D*lanXc^2)+(GG*lanYb^2+DD*lanXb^2)+a8*lanY5^2 -B*lanXc -(GG-DD)*lanXb*lanYb -(GG*lanYb^2+DD*lanXb^2) 0;
           C 0 B*lanYc -B*lanXc A+AA BB*lanYb -BB*lanXb CC;
           0 B*lanYc -(GG*lanXb^2+DD*lanYb^2) -(GG-DD)*lanXb*lanYb BB*lanYb (G*lanXc^2+D*lanYc^2)+(GG*lanXb^2+DD*lanYb^2)+a5*lanX4^2 (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb+a6*lanX4*lanY4 B*lanYc;
           0 -B*lanXc -(GG-DD)*lanXb*lanYb -(GG*lanYb^2+DD*lanXb^2) -BB*lanXb (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb+a5*lanX4*lanY4 (G*lanYc^2+D*lanXc^2)+(GG*lanYb^2+DD*lanXb^2)+a6*lanY4^2 -B*lanXc;
           0 C 0 0 CC B*lanYc -B*lanXc A+AA];
    % Define the applied load
F = [P3;P6;P7*i;P8;P9;P10;P11;P12];
it = 0; % initialize iteration count
residual = 100; % initialize residual
% calculate Force
while (residual > tolerance)
% assemble global K matrix
a1=((Rkic-Rkpc)/((1+(abs((Rkic*u(1))/Myc))^nc)^(1/nc)))+Rkpc;
a2=((Rkic-Rkpc)/((1+(abs((Rkic*u(2))/Myc))^nc)^(1/nc)))+Rkpc;
a5=((Rki-Rkp)/((1+(abs((Rki*u(6))/Fyb))^n)^(1/n)))+Rkp;
a6=((Rki-Rkp)/((1+(abs((Rki*u(7))/Fyb))^n)^(1/n)))+Rkp;
a7=((Rki-Rkp)/((1+(abs((Rki*u(3))/Fyb))^n)^(1/n)))+Rkp;
a8=((Rki-Rkp)/((1+(abs((Rki*u(4))/Fyb))^n)^(1/n)))+Rkp;
        %  3    6   7    8    9    10     11    12
      K = [A+a1 0 B*lanYc -B*lanXc C 0 0 0;
           0 A+a2 0 0 0 B*lanYc -B*lanXc C;
           B*lanYc 0 (G*lanXc^2+D*lanYc^2)+(GG*lanXb^2+DD*lanYb^2)+a7*lanX5^2 (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb B*lanYc+a8*lanX5*lanY5 -(GG*lanXb^2+DD*lanYb^2) -(GG-DD)*lanXb*lanYb 0;
          -B*lanXc 0 (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb+a7*lanX5*lanY5 (G*lanYc^2+D*lanXc^2)+(GG*lanYb^2+DD*lanXb^2)+a8*lanY5^2 -B*lanXc -(GG-DD)*lanXb*lanYb -(GG*lanYb^2+DD*lanXb^2) 0;
           C 0 B*lanYc -B*lanXc A+AA BB*lanYb -BB*lanXb CC;
           0 B*lanYc -(GG*lanXb^2+DD*lanYb^2) -(GG-DD)*lanXb*lanYb BB*lanYb (G*lanXc^2+D*lanYc^2)+(GG*lanXb^2+DD*lanYb^2)+a5*lanX4^2 (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb+a6*lanX4*lanY4 B*lanYc;
           0 -B*lanXc -(GG-DD)*lanXb*lanYb -(GG*lanYb^2+DD*lanXb^2) -BB*lanXb (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb+a5*lanX4*lanY4 (G*lanYc^2+D*lanXc^2)+(GG*lanYb^2+DD*lanXb^2)+a6*lanY4^2 -B*lanXc;
           0 C 0 0 CC B*lanYc -B*lanXc A+AA];
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
    ucb3(i)=u(3);Hcb3(i)=(((Rki-Rkp)*ucb3(i))/((1+(abs((Rki*ucb3(i))/Fyb))^n)^(1/n)))+Rkp*ucb3(i);
   if it < itermax
   fprintf('(+)It is converged in %1.0f iterations for increment %1.0f\n',it,i)
   end
      if 1 == isnan(u) ;break;end % Check for Softening
 if abs(u(1)) >= tuc;disp('      ## spring at support reached to Ultimate Rotation ##');break;end
 if or(abs(u(3))>= dub,abs(u(6))>= dub);disp('      ## Brace displacement reached to Ultimate displacement ##');break;end
end
D1=[0;U7'];F1=[0;F7'];
% Structural Lateral Stiffness [DOF(7)]
s=size(D1,1);for i=1:s-1;Kz1(i)=(F1(i+1)-F1(i))/(D1(i+1)-D1(i));end
% Coordinate of element
corx=[0 0 L L];cory=[0 h h 0];
% For element - Brace
corx4=[0;L];cory4=[0;h];corx5=[L;0];cory5=[0;h];
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
% For element 4 - Brace
X4N1=[0;L+u(6)];Y4N1=[0;h+u(7)];
% For element 5 - Brace
X5N1=[L;0+u(3)];Y5N1=[0;h+u(4)];
%% Second order Nonlinear Analysis
disp('#################################################');
disp('#        Second-order Nonlinear Analysis        #');
disp('#################################################');
U = zeros(8,1);% initial guess value
% Element stiffness coefficient for element 1 and 2
A=(4*EIc/h)+(2*P8*h/15);B=(6*EIc/h^2)+(P8/10);C=(2*EIc/h)-(P8*h/30);D=(12*EIc/h^3)+(6*P8/(5*h));G=(EAc/h)+(P8/h);
% Element stiffness coefficient for element 3
AA=(4*EIb/L);BB=(6*EIb/L^2);CC=(2*EIb/L);DD=(12*EIb/h^3);GG=(EAb/L);
for i=1:m
a1=((Rkic-Rkpc)/((1+(abs((Rkic*U(1))/Myc))^nc)^(1/nc)))+Rkpc;
a2=((Rkic-Rkpc)/((1+(abs((Rkic*U(2))/Myc))^nc)^(1/nc)))+Rkpc;
a5=((Rki-Rkp)/((1+(abs((Rki*U(6))/Fyb))^n)^(1/n)))+Rkp;
a6=((Rki-Rkp)/((1+(abs((Rki*U(7))/Fyb))^n)^(1/n)))+Rkp;
a7=((Rki-Rkp)/((1+(abs((Rki*U(3))/Fyb))^n)^(1/n)))+Rkp;
a8=((Rki-Rkp)/((1+(abs((Rki*U(4))/Fyb))^n)^(1/n)))+Rkp;
% initial assemble global K matrix
        %  3    6   7    8    9    10     11    12
  Kinit = [A+a1 0 B*lanYc -B*lanXc C 0 0 0;
           0 A+a2 0 0 0 B*lanYc -B*lanXc C;
           B*lanYc 0 (G*lanXc^2+D*lanYc^2)+(GG*lanXb^2+DD*lanYb^2)+a7*lanX5^2 (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb B*lanYc+a8*lanX5*lanY5 -(GG*lanXb^2+DD*lanYb^2) -(GG-DD)*lanXb*lanYb 0;
          -B*lanXc 0 (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb+a7*lanX5*lanY5 (G*lanYc^2+D*lanXc^2)+(GG*lanYb^2+DD*lanXb^2)+a8*lanY5^2 -B*lanXc -(GG-DD)*lanXb*lanYb -(GG*lanYb^2+DD*lanXb^2) 0;
           C 0 B*lanYc -B*lanXc A+AA BB*lanYb -BB*lanXb CC;
           0 B*lanYc -(GG*lanXb^2+DD*lanYb^2) -(GG-DD)*lanXb*lanYb BB*lanYb (G*lanXc^2+D*lanYc^2)+(GG*lanXb^2+DD*lanYb^2)+a5*lanX4^2 (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb+a6*lanX4*lanY4 B*lanYc;
           0 -B*lanXc -(GG-DD)*lanXb*lanYb -(GG*lanYb^2+DD*lanXb^2) -BB*lanXb (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb+a5*lanX4*lanY4 (G*lanYc^2+D*lanXc^2)+(GG*lanYb^2+DD*lanXb^2)+a6*lanY4^2 -B*lanXc;
           0 C 0 0 CC B*lanYc -B*lanXc A+AA];
% Define the applied load
F = [P3;P6;P7*i;P8;P9;P10;P11;P12];
it = 0; % initialize iteration count
residual = 100; % initialize residual
% calculate Force
while (residual > tolerance)
% assemble global K matrix
a1=((Rkic-Rkpc)/((1+(abs((Rkic*U(1))/Myc))^nc)^(1/nc)))+Rkpc;
a2=((Rkic-Rkpc)/((1+(abs((Rkic*U(2))/Myc))^nc)^(1/nc)))+Rkpc;
a5=((Rki-Rkp)/((1+(abs((Rki*U(6))/Fyb))^n)^(1/n)))+Rkp;
a6=((Rki-Rkp)/((1+(abs((Rki*U(7))/Fyb))^n)^(1/n)))+Rkp;
a7=((Rki-Rkp)/((1+(abs((Rki*U(3))/Fyb))^n)^(1/n)))+Rkp;
a8=((Rki-Rkp)/((1+(abs((Rki*U(4))/Fyb))^n)^(1/n)))+Rkp;
        %  3    6   7    8    9    10     11    12
      K = [A+a1 0 B*lanYc -B*lanXc C 0 0 0;
           0 A+a2 0 0 0 B*lanYc -B*lanXc C;
           B*lanYc 0 (G*lanXc^2+D*lanYc^2)+(GG*lanXb^2+DD*lanYb^2)+a7*lanX5^2 (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb B*lanYc+a8*lanX5*lanY5 -(GG*lanXb^2+DD*lanYb^2) -(GG-DD)*lanXb*lanYb 0;
          -B*lanXc 0 (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb+a7*lanX5*lanY5 (G*lanYc^2+D*lanXc^2)+(GG*lanYb^2+DD*lanXb^2)+a8*lanY5^2 -B*lanXc -(GG-DD)*lanXb*lanYb -(GG*lanYb^2+DD*lanXb^2) 0;
           C 0 B*lanYc -B*lanXc A+AA BB*lanYb -BB*lanXb CC;
           0 B*lanYc -(GG*lanXb^2+DD*lanYb^2) -(GG-DD)*lanXb*lanYb BB*lanYb (G*lanXc^2+D*lanYc^2)+(GG*lanXb^2+DD*lanYb^2)+a5*lanX4^2 (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb+a6*lanX4*lanY4 B*lanYc;
           0 -B*lanXc -(GG-DD)*lanXb*lanYb -(GG*lanYb^2+DD*lanXb^2) -BB*lanXb (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb+a5*lanX4*lanY4 (G*lanYc^2+D*lanXc^2)+(GG*lanYb^2+DD*lanXb^2)+a6*lanY4^2 -B*lanXc;
           0 C 0 0 CC B*lanYc -B*lanXc A+AA];
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
    Ucb3(i)=U(3);HCB3(i)=(((Rki-Rkp)*Ucb3(i))/((1+(abs((Rki*Ucb3(i))/Fyb))^n)^(1/n)))+Rkp*Ucb3(i);
   if it < itermax
   fprintf('(+)It is converged in %1.0f iterations for increment %1.0f\n',it,i)
   end
   if 1 == isnan(U) ;break;end % Check for Softening
   if abs(U(1)) >= tuc;disp('      ## spring at support reached to Ultimate Rotation ##');break;end
   if or(abs(U(3))>= dub,abs(U(6))>= dub);disp('      ## Brace displacement reached to Ultimate displacement ##');break;end
end
D2=[0;UI7'];F2=[0;FI7'];
% Structural Lateral Stiffness [DOF(7)]
s=size(D2,1);for i=1:s-1;Kz2(i)=(F2(i+1)-F2(i))/(D2(i+1)-D2(i));end
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
% For element 4 - Brace
X4N2=[0;L+U(6)];Y4N2=[0;h+U(7)];
% For element 5 - Brace
X5N2=[L;0+U(3)];Y5N2=[0;h+U(4)];
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
for i=1:200;
    Uc=tuc/200;uc(i)=Uc*i;
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
%% Brace - Hinge Ploting - Displacement
for i=1:30;
    Ucbi=dub/30;ucbi(i)=Ucbi*i;
    Hcbi(i)=(((Rki-Rkp)*ucbi(i))/((1+(abs((Rki*ucbi(i))/Fyb))^n)^(1/n)))+Rkp*ucbi(i);
end
ucbi =[0 ucbi];Hcbi =[0 Hcbi];
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
fprintf('=== 1st-Order Nonlinear ==+== 2nd-Order Nonlinear =====\n');
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
IMAGE=imread('PushoverAnalOfNonlinSemiRigidSpringWithElaFrameUnderDistriLoadSTIFFNESS3.jpg');
image(IMAGE);axis image;axis off;
 title('Stiffness Matrix for element 4 & 5 - Brace','FontSize',12,'color','b')
figure(4)
IMAGE=imread('PushoverAnalOfBraceTwo.jpg');
image(IMAGE);axis image;axis off;
figure(5)
IMAGE=imread('PushoverAnalOfBraceNonlinSemiRigSpriUnderDistriLoadl.jpg');
image(IMAGE);axis image;axis off;
figure(6)
IMAGE=imread('Modified Newton Raphson Method for single increment.jpg');
image(IMAGE);axis image;axis off;
figure(7)
p1=plot(I1,DU1,'black',I2,DU2,'--blue');grid on;set(p1,'LineWidth',2);
xlabel('increment');ylabel('Residual');
title('Residual-increment diagram','color','b');
legend('1st-Order Nonlinear Analysis','2nd-Order Nonlinear Analysis','Location','NorthEastOutside');
figure(8)
p1=plot(I1,IT1,'black',I2,IT2,'--blue');grid on;set(p1,'LineWidth',2);
xlabel('increment');ylabel('Iteration');
title('Iteration-increment diagram','color','b');
legend('1st-Order Nonlinear Analysis','2nd-Order Nonlinear Analysis','Location','NorthEastOutside');
figure(9)
p1=plot(uc,Hc,'b',Xuc,YHc,'--r');grid on;set(p1,'LineWidth',3);
legend('Column Hinge','Column Hinge Curve bilinear fitted','Location','NorthEastOutside');
xlabel('Rotation (rad)');ylabel('Moment (kN.mm)');
title('Moment-Rotation Relation of Nonlin. Semi-Rigid Connection','color','b')
figure(10)
p1=plot(ucbi,Hcbi,'b');grid on;set(p1,'LineWidth',3);
xlabel('Displacement (mm)');ylabel('Force (kN)');
title('Nonlin. Behavior of  Brace Force-Displacement Relation','color','b')
figure(11)
p1=plot([0 uc1],[0 Hc1],[0 uc2],[0 Hc2],[0 Uc1],[0 HC1],'--',[0 Uc2],[0 HC2],'--');grid on;set(p1,'LineWidth',3);
legend('1st-Order anal. Hing DOF(3)','1st-Order anal. Hing DOF(6)','2nd-Order anal. Hing DOF(3)','2nd-Order anal. Hing DOF(6)','Location','NorthEastOutside');
xlabel('Rotation (rad)');ylabel('Moment (kN.mm)');
title('Moment-Rotation Diagram of Connection(Spring) During the Analysis','color','b');
figure(12)
p1=plot([0 ucb3],[0 Hcb3],'--b',[0 Ucb3],[0 HCB3],'-.r');grid on;set(p1,'LineWidth',3);
legend('1st-Order Brace behavior','2nd-Order Brace behavior','Location','NorthEastOutside');
xlabel('Displacement (mm)');ylabel('Force (kN)');
title('Nonlin. behavior of  Brace Force-Displacement Relation During the Analysis','color','b')
figure(13)
p1=plot(D1,F1,'b',X1,Y1,'--r',D2,F2,'g',X2,Y2,'--');grid on;set(p1,'LineWidth',3);
legend('1st-Order Nonlinear Analysis','1st-Order Nonlinear Curve fitted','2nd-Order Nonlinear Analysis','2nd-Order Nonlinear Curve fitted','Location','NorthEastOutside');
xlabel('Displacement (mm) [DOF (7)]');ylabel('Force (kN) [DOF (7)]');
title(['Force-Displacement Diagram of Frame with Nonlin. Semi-Rigid Connection',' 1-order Nonlin. Ductility Rito: ',num2str(Z3),'   2-order Nonlin. Ductility Rito: ',num2str(Z4)],'color','b')
figure(14)
p1=plot(U7,Kz1,'black',UI7,Kz2,'blue--');grid on;set(p1,'LineWidth',3);
legend('1st-Order Nonlinear Analysis','2nd-Order Nonlinear Analysis','Location','NorthEastOutside');
xlabel('Displacement (mm) [DOF (7)]');ylabel('Structural Lateral Stiffness (kN/mm) [DOF (7)]');
title('Structural Lateral Stiffness-Displacement Diagram of Elastic Frame with Nonlin. Semi-Rigid Connection','color','b')
figure(15)
p1=plot(corx,cory,'--b',corx4,cory4,'--b',corx5,cory5,'--b',X1N1,Y1N1,'--r',X2N1,Y2N1,'--r',X3N1,Y3N1,'--r',X4N1,Y4N1,'--r',X5N1,Y5N1,'--r',X1N2,Y1N2,'-.g',X2N2,Y2N2,'-.g',X3N2,Y3N2,'-.g',X4N2,Y4N2,'-.g',X5N2,Y5N2,'-.g');grid on;set(p1,'LineWidth',3);
legend('BLUE: Not Loading','RED: 1st-Order Nonlinear Analysis','GREEN: 2nd-Order Nonlinear Analysis','Location','NorthOutside');
xlabel('X (mm)');ylabel('Y (mm)');
title('Small Deflection Theory - Last Step Displacement of Nonlin. behavior of Brace and Elastic Frame with Nonlin. Semi-Rigid Connection','color','b')
%%  print time of computation
totaltime = cputime - starttime;
fprintf('\nTotal time (s)= %7.4f \n\n',totaltime)
