%***********************************************************%
%                >> IN THE NAME OF GOD <<                   %
% Analysis of 1st order and 2nd order Nonlinear Semi Rigid  %
% Connection Frame subjected to Pushover lateral load       %
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
EIc = 200*100^4/12; % [kN.mm^2] column
EAc = 200*10000; % [kN]
EIb = 200*50^4/12; % [kN.mm^2] beam
EAb = 200*(50)^2; % [kN]
m = 1000; % number of calculation (Load Steps)
itermax = 4000;% maximum number of iterations
tolerance = 1e-9; % specified tolerance for convergence
u = zeros(8,1);% initial guess value
P3 =0; % [kN.mm] Moment [DOF (3)]
P6 =0; % [kN.mm] Moment [DOF (6)]
P7=0.1; % [kN] Horizental Force [DOF (7)] Incremantal Loading
P8=-.5*W*L; % [kN] Vertical Force [DOF (8)]
P9 =-(W*L^2)/12; % [kN.mm] Moment [DOF (9)]
P10=0; % [kN] Horizental Force [DOF (10)]
P11=-.5*W*L; % [kN] Vertical Force [DOF (11)]
P12=+(W*L^2)/12; % [kN.mm] Moment [DOF (12)]
lanXc=0;lanYc=1;
lanXb=1;lanYb=0;
%% Nonlinear Rotational Spring of columns
tyc=.001; % Yield rotaion
Myc=40e+3; % Yield moment
tuc=.025;  % Ultimate rotation
Muc=1.25*Myc;  % Ultimate moment
nc = 9; % Moment-rotation shape parameter
Rkic=Myc/tyc;
Rkpc=(Muc-Myc)/(tuc-tyc);
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
a1=((Rkic-Rkpc)/((1+(abs((Rkic*u(1))/Myc))^nc)^(1/nc)))+Rkpc;
a2=((Rkic-Rkpc)/((1+(abs((Rkic*u(2))/Myc))^nc)^(1/nc)))+Rkpc;
% initial assemble global K matrix
        %  3    6   7    8    9    10     11    12
  Kinit = [A+a1 0 B*lanYc -B*lanXc C 0 0 0;
           0 A+a2 0 0 0 B*lanYc -B*lanXc C;
           B*lanYc 0 (G*lanXc^2+D*lanYc^2)+(GG*lanXb^2+DD*lanYb^2) (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb B*lanYc -(GG*lanXb^2+DD*lanYb^2) -(GG-DD)*lanXb*lanYb 0;
          -B*lanXc 0 (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb (G*lanYc^2+D*lanXc^2)+(GG*lanYb^2+DD*lanXb^2) -B*lanXc -(GG-DD)*lanXb*lanYb -(GG*lanYb^2+DD*lanXb^2) 0;
           C 0 B*lanYc -B*lanXc A+AA BB*lanYb -BB*lanXb CC;
           0 B*lanYc -(GG*lanXb^2+DD*lanYb^2) -(GG-DD)*lanXb*lanYb BB*lanYb (G*lanXc^2+D*lanYc^2)+(GG*lanXb^2+DD*lanYb^2) (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb B*lanYc;
           0 -B*lanXc -(GG-DD)*lanXb*lanYb -(GG*lanYb^2+DD*lanXb^2) -BB*lanXb (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb (G*lanYc^2+D*lanXc^2)+(GG*lanYb^2+DD*lanXb^2) -B*lanXc;
           0 C 0 0 CC B*lanYc -B*lanXc A+AA];
% Define the applied load
F = [P3;P6;P7*i;P8;P9;P10;P11;P12];
it = 0; % initialize iteration count
residual = 100; % initialize residual
% calculate Force
while (residual > tolerance)
% assemble global K matrix
a1=((Rkic-Rkpc)/((1+(abs((Rkic*u(1))/Myc))^nc)^(1/nc)))+Rkpc;% Column hinge [DOF (3)]
a2=((Rkic-Rkpc)/((1+(abs((Rkic*u(2))/Myc))^nc)^(1/nc)))+Rkpc;% Column hinge [DOF (6)]
        %  3    6   7    8    9    10     11    12
      K = [A+a1 0 B*lanYc -B*lanXc C 0 0 0;
           0 A+a2 0 0 0 B*lanYc -B*lanXc C;
           B*lanYc 0 (G*lanXc^2+D*lanYc^2)+(GG*lanXb^2+DD*lanYb^2) (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb B*lanYc -(GG*lanXb^2+DD*lanYb^2) -(GG-DD)*lanXb*lanYb 0;
          -B*lanXc 0 (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb (G*lanYc^2+D*lanXc^2)+(GG*lanYb^2+DD*lanXb^2) -B*lanXc -(GG-DD)*lanXb*lanYb -(GG*lanYb^2+DD*lanXb^2) 0;
           C 0 B*lanYc -B*lanXc A+AA BB*lanYb -BB*lanXb CC;
           0 B*lanYc -(GG*lanXb^2+DD*lanYb^2) -(GG-DD)*lanXb*lanYb BB*lanYb (G*lanXc^2+D*lanYc^2)+(GG*lanXb^2+DD*lanYb^2) (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb B*lanYc;
           0 -B*lanXc -(GG-DD)*lanXb*lanYb -(GG*lanYb^2+DD*lanXb^2) -BB*lanXb (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb (G*lanYc^2+D*lanXc^2)+(GG*lanYb^2+DD*lanXb^2) -B*lanXc;
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
 % Internal element force           
         % Displacement Transformation Matrix
        T = [lanXc lanYc 0 0 0 0;
            -lanYc lanXc 0 0 0 0;
                   0 0 1 0 0 0;
             0 0 0 lanXc lanYc 0;
            0 0 0 -lanYc lanXc 0;
                   0 0 0 0 0 1];
        % Stiffness Matrix for each element
        Kele1 = [G 0 0 -G 0 0;
                 0 D B 0 -D B;
              0 B A 0 -B C;
                 -G 0 0 G 0 0;
               0 -D -B 0 D -B;
                 0 B C 0 -B A]; 
        Fele1 = Kele1*T*[0;0;u(1);u(3);u(4);u(5)];% Internal Force of element column 1
        Kele2 = [G 0 0 -G 0 0;
                 0 D B 0 -D B;
              0 B A 0 -B C;
                 -G 0 0 G 0 0;
               0 -D -B 0 D -B;
                 0 B C 0 -B A];        
        Fele2 = Kele2*T*[0;0;u(2);u(6);u(7);u(8)];% Internal Force of element column 2
       T3 = [lanXb lanYb 0 0 0 0;
            -lanYb lanXb 0 0 0 0;
                     0 0 1 0 0 0;
             0 0 0 lanXb lanYb 0;
            0 0 0 -lanYb lanXb 0;
                     0 0 0 0 0 1];
        % Stiffness Matrix for each element
        Kele3 = [GG 0 0 -GG 0 0;
               0 DD BB 0 -DD BB;
               0 BB AA 0 -BB CC;
                 -GG 0 0 GG 0 0;
             0 -DD -BB 0 DD -BB;
               0 BB CC 0 -BB AA]; 
        Fele3 = Kele3*T3*[u(3);u(4);u(5);u(6);u(7);u(8)]+[0;.5*W*L;+(W*L^2)/12;0;.5*W*L;-(W*L^2)/12];% Internal Force of element beam
% Internal Force and Dispalcement for each increment
    INT1_N1_f1(i) = roundn(Fele1(1),-3);INT2_N1_f1(i) = roundn(Fele2(1),-3);INT3_N1_f1(i) = roundn(Fele3(1),-3);
    INT1_N1_f2(i) = roundn(Fele1(2),-3);INT2_N1_f2(i) = roundn(Fele2(2),-3);INT3_N1_f2(i) = roundn(Fele3(2),-3);
    INT1_N1_f3(i) = roundn(Fele1(3),-3);INT2_N1_f3(i) = roundn(Fele2(3),-3);INT3_N1_f3(i) = roundn(Fele3(3),-3);
    INT1_N1_f4(i) = roundn(Fele1(4),-3);INT2_N1_f4(i) = roundn(Fele2(4),-3);INT3_N1_f4(i) = roundn(Fele3(4),-3);
    INT1_N1_f5(i) = roundn(Fele1(5),-3);INT2_N1_f5(i) = roundn(Fele2(5),-3);INT3_N1_f5(i) = roundn(Fele3(5),-3);
    INT1_N1_f6(i) = roundn(Fele1(6),-3);INT2_N1_f6(i) = roundn(Fele2(6),-3);INT3_N1_f6(i) = roundn(Fele3(6),-3);    
    TBSN1(i)=[INT1_N1_f2(i)+INT2_N1_f2(i)]; % Total Base Shear of Columns
    TMSN1(i)=[INT1_N1_f3(i)+INT2_N1_f3(i)]; % Total Base Moment of Columns
    F7(i) = F(3);
    U3(i) = u(1);U6(i) = u(2);
    U7(i) = u(3);U8(i) = u(4);
    U9(i) = u(5);U10(i) = u(6);
    U11(i) = u(7);U12(i) = u(8);
    DU1(i)=residual;I1(i)=i;IT1(i)=it;
    % Moment-Rotation of spring during the analysis
    uc1(i)=u(1);Hc1(i)=(((Rkic-Rkpc)*uc1(i))/((1+(abs((Rkic*uc1(i))/Myc))^nc)^(1/nc)))+Rkpc*uc1(i);
    uc2(i)=u(2);Hc2(i)=(((Rkic-Rkpc)*uc2(i))/((1+(abs((Rkic*uc2(i))/Myc))^nc)^(1/nc)))+Rkpc*uc2(i);
   if it < itermax
   fprintf('(+)It is converged in %1.0f iterations for increment %1.0f\n',it,i)
   end
      if 1 == isnan(u) ;break;end % Check for Softening
      if abs(u(1)) >= tuc;disp('      ## spring at support [DOF (3)] reached to Ultimate Rotation ##');break;end
      if abs(u(2)) >= tuc;disp('      ## spring at support [DOF (6)] reached to Ultimate Rotation ##');break;end
end
D1=[0;U7'];F1=[0;F7'];% Force-displacement of 1st-order Nonlinear analysis
% Structural Lateral Stiffness [DOF(7)]
s=size(D1,1);for i=1:s-1;Kz1(i)=(F1(i+1)-F1(i))/(D1(i+1)-D1(i));end
% Coordinate of element
corx=[0 0 L L];cory=[0 h h 0];
% displaying shape of deflection along the length for First-order Nonlinear Analysis
% For element 1
for j=1:21;
    y1(j)=((h+u(4))/20)*(j-1);
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
for j=1:21;
    y3(j)=((h+u(7))/20)*(j-1);
    N1 =(1/h^3)*(h^3-3*h*y3(j)^2+2*y3(j)^3);
    N2 =(1/h^2)*(y3(j)*h^2-2*h*y3(j)^2+y3(j)^3);
    N3 =(1/h^3)*(3*h*y3(j)^2-2*y3(j)^3);
    N4 =(1/h^2)*(-h*y3(j)^2+y3(j)^3);
    x3(j) =[N1 N2 N3 N4]*[0;u(2);u(6);u(8)]+L;
end
X1N1=x1;Y1N1=y1;X2N1=x2;Y2N1=y2;X3N1=x3;Y3N1=y3;
%% First-order Linear Analysis
for i=1:i 
% calculate Force
 F = [P7*i;P9;P12];
% assemble global K matrix
        K=[D+D B B;
           B A+AA C;
           B C A+AA];
        uL1 = K^-1 *F;       
% Force and Dispalcement for each increment
    F1L1(i) = F(1);
    u1L(i) = uL1(1);
    u2L(i) = uL1(2);
    u3L(i) = uL1(3);
end
D1L=[0;u1L'];F1L=[0;F1L1'];% Force-displacement of 1st-order Linear analysis
%% Second order Nonlinear Analysis
disp('#################################################');
disp('#        Second-order Nonlinear Analysis        #');
disp('#################################################');
U = zeros(8,1);% initial guess value
% Element stiffness coefficient for element 1 and 2
A=(4*EIc/h)+(2*P8*h/15);B=(6*EIc/h^2)+(P8/10);C=(2*EIc/h)-(P8*h/30);D=(12*EIc/h^3)+(6*P8/(5*h));G=(EAc/h)+(P8/h);
% Element stiffness coefficient for element 3
AA=(4*EIb/L);BB=(6*EIb/L^2);CC=(2*EIb/L);DD=(12*EIb/h^3);GG=(EAb/L);
for ii=1:m
a1=((Rkic-Rkpc)/((1+(abs((Rkic*U(1))/Myc))^nc)^(1/nc)))+Rkpc;
a2=((Rkic-Rkpc)/((1+(abs((Rkic*U(2))/Myc))^nc)^(1/nc)))+Rkpc;
% initial assemble global K matrix
        %  3    6   7    8    9    10     11    12
  Kinit = [A+a1 0 B*lanYc -B*lanXc C 0 0 0;
           0 A+a2 0 0 0 B*lanYc -B*lanXc C;
           B*lanYc 0 (G*lanXc^2+D*lanYc^2)+(GG*lanXb^2+DD*lanYb^2) (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb B*lanYc -(GG*lanXb^2+DD*lanYb^2) -(GG-DD)*lanXb*lanYb 0;
          -B*lanXc 0 (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb (G*lanYc^2+D*lanXc^2)+(GG*lanYb^2+DD*lanXb^2) -B*lanXc -(GG-DD)*lanXb*lanYb -(GG*lanYb^2+DD*lanXb^2) 0;
           C 0 B*lanYc -B*lanXc A+AA BB*lanYb -BB*lanXb CC;
           0 B*lanYc -(GG*lanXb^2+DD*lanYb^2) -(GG-DD)*lanXb*lanYb BB*lanYb (G*lanXc^2+D*lanYc^2)+(GG*lanXb^2+DD*lanYb^2) (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb B*lanYc;
           0 -B*lanXc -(GG-DD)*lanXb*lanYb -(GG*lanYb^2+DD*lanXb^2) -BB*lanXb (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb (G*lanYc^2+D*lanXc^2)+(GG*lanYb^2+DD*lanXb^2) -B*lanXc;
           0 C 0 0 CC B*lanYc -B*lanXc A+AA];
% Define the applied load
F = [P3;P6;P7*ii;P8;P9;P10;P11;P12];
it = 0; % initialize iteration count
residual = 100; % initialize residual
% calculate Force
while (residual > tolerance)
% assemble global K matrix
a1=((Rkic-Rkpc)/((1+(abs((Rkic*U(1))/Myc))^nc)^(1/nc)))+Rkpc;
a2=((Rkic-Rkpc)/((1+(abs((Rkic*U(2))/Myc))^nc)^(1/nc)))+Rkpc;
        %  3    6   7    8    9    10     11    12
      K = [A+a1 0 B*lanYc -B*lanXc C 0 0 0;
           0 A+a2 0 0 0 B*lanYc -B*lanXc C;
           B*lanYc 0 (G*lanXc^2+D*lanYc^2)+(GG*lanXb^2+DD*lanYb^2) (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb B*lanYc -(GG*lanXb^2+DD*lanYb^2) -(GG-DD)*lanXb*lanYb 0;
          -B*lanXc 0 (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb (G*lanYc^2+D*lanXc^2)+(GG*lanYb^2+DD*lanXb^2) -B*lanXc -(GG-DD)*lanXb*lanYb -(GG*lanYb^2+DD*lanXb^2) 0;
           C 0 B*lanYc -B*lanXc A+AA BB*lanYb -BB*lanXb CC;
           0 B*lanYc -(GG*lanXb^2+DD*lanYb^2) -(GG-DD)*lanXb*lanYb BB*lanYb (G*lanXc^2+D*lanYc^2)+(GG*lanXb^2+DD*lanYb^2) (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb B*lanYc;
           0 -B*lanXc -(GG-DD)*lanXb*lanYb -(GG*lanYb^2+DD*lanXb^2) -BB*lanXb (G-D)*lanXc*lanYc+(GG-DD)*lanXb*lanYb (G*lanYc^2+D*lanXc^2)+(GG*lanYb^2+DD*lanXb^2) -B*lanXc;
           0 C 0 0 CC B*lanYc -B*lanXc A+AA];
        f=F-K*U;
        %calculate du
        du = Kinit^-1 *(f);
        U = U+du; % update u
        %Calculate the residual (internal-external force deviation)
        residual = max(abs(du)); % evaluate residual
        it = it + 1; % increment iteration count
        if it == itermax
          fprintf('(-)For increment %1.0f trail iteration reached to Ultimate %1.0f\n',ii,it)
             disp('    ## The solution for this step is not converged ##')
        end
end
%% Force and Dispalcement for each increment
% Internal element force           
         % Displacement Transformation Matrix
        Tc = [lanXc lanYc 0 0 0 0;
            -lanYc lanXc 0 0 0 0;
                   0 0 1 0 0 0;
             0 0 0 lanXc lanYc 0;
            0 0 0 -lanYc lanXc 0;
                   0 0 0 0 0 1];
        % Stiffness Matrix for each element
        Kele1 = [G 0 0 -G 0 0;
                 0 D B 0 -D B;
              0 B A 0 -B C;
                 -G 0 0 G 0 0;
               0 -D -B 0 D -B;
                 0 B C 0 -B A]; 
        Fele1 = Kele1*Tc*[0;0;U(1);U(3);U(4);U(5)];% Internal Force of element column 1
        Kele2 = [G 0 0 -G 0 0;
                 0 D B 0 -D B;
              0 B A 0 -B C;
                 -G 0 0 G 0 0;
               0 -D -B 0 D -B;
                 0 B C 0 -B A];        
        Fele2 = Kele2*Tc*[0;0;U(2);U(6);U(7);U(8)];% Internal Force of element column 2
       T3 = [lanXb lanYb 0 0 0 0;
            -lanYb lanXb 0 0 0 0;
                     0 0 1 0 0 0;
             0 0 0 lanXb lanYb 0;
            0 0 0 -lanYb lanXb 0;
                     0 0 0 0 0 1];
        % Stiffness Matrix for each element
        Kele3 = [GG 0 0 -GG 0 0;
               0 DD BB 0 -DD BB;
               0 BB AA 0 -BB CC;
                 -GG 0 0 GG 0 0;
             0 -DD -BB 0 DD -BB;
               0 BB CC 0 -BB AA]; 
        Fele3 = Kele3*T3*[U(3);U(4);U(5);U(6);U(7);U(8)]+[0;.5*W*L;+(W*L^2)/12;0;.5*W*L;-(W*L^2)/12];% Internal Force of element beam
% Internal Force and Dispalcement for each increment
    INT1_N2_f1(ii) = roundn(Fele1(1),-3);INT2_N2_f1(ii) = roundn(Fele2(1),-3);INT3_N2_f1(ii) = roundn(Fele3(1),-3);
    INT1_N2_f2(ii) = roundn(Fele1(2),-3);INT2_N2_f2(ii) = roundn(Fele2(2),-3);INT3_N2_f2(ii) = roundn(Fele3(2),-3);
    INT1_N2_f3(ii) = roundn(Fele1(3),-3);INT2_N2_f3(ii) = roundn(Fele2(3),-3);INT3_N2_f3(ii) = roundn(Fele3(3),-3);
    INT1_N2_f4(ii) = roundn(Fele1(4),-3);INT2_N2_f4(ii) = roundn(Fele2(4),-3);INT3_N2_f4(ii) = roundn(Fele3(4),-3);
    INT1_N2_f5(ii) = roundn(Fele1(5),-3);INT2_N2_f5(ii) = roundn(Fele2(5),-3);INT3_N2_f5(ii) = roundn(Fele3(5),-3);
    INT1_N2_f6(ii) = roundn(Fele1(6),-3);INT2_N2_f6(ii) = roundn(Fele2(6),-3);INT3_N2_f6(ii) = roundn(Fele3(6),-3);    
    TBSN2(ii)=[INT1_N2_f2(ii)+INT2_N2_f2(ii)]; % Total Base Shear of Columns
    TMSN2(ii)=[INT1_N2_f3(ii)+INT2_N2_f3(ii)]; % Total Base Moment of Columns
    FI7(ii) = F(3);
    UI3(ii) = U(1);UI6(ii) = U(2);
    UI7(ii) = U(3);UI8(ii) = U(4);
    UI9(ii) = U(5);UI10(ii) = U(6);
    UI11(ii) = U(7);UI12(ii) = U(8);
    DU2(ii)=residual;I2(ii)=ii;IT2(ii)=it;
    % Moment-Rotation of spring during the analysis
    Uc1(ii)=U(1);HC1(ii)=(((Rkic-Rkpc)*Uc1(ii))/((1+(abs((Rkic*Uc1(ii))/Myc))^nc)^(1/nc)))+Rkpc*Uc1(ii);
    Uc2(ii)=U(2);HC2(ii)=(((Rkic-Rkpc)*Uc2(ii))/((1+(abs((Rkic*Uc2(ii))/Myc))^nc)^(1/nc)))+Rkpc*Uc2(ii);
   if it < itermax
   fprintf('(+)It is converged in %1.0f iterations for increment %1.0f\n',it,ii)
   end
   if 1 == isnan(U) ;break;end % Check for Softening
      if abs(U(1)) >= tuc;disp('      ## spring at support [DOF (3)] reached to Ultimate Rotation ##');break;end
      if abs(U(2)) >= tuc;disp('      ## spring at support [DOF (6)] reached to Ultimate Rotation ##');break;end
end
D2=[0;UI7'];F2=[0;FI7'];% Force-displacement of 2nd-order Nonlinear analysis
% Structural Lateral Stiffness [DOF(7)]
s=size(D2,1);for i=1:s-1;Kz2(i)=(F2(i+1)-F2(i))/(D2(i+1)-D2(i));end
% displaying shape of deflection along the length for First-order Nonlinear Analysis
% For element 1
for j=1:21;
    y1(j)=((h+U(4))/20)*(j-1);
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
for j=1:21;
    y3(j)=((h+U(7))/20)*(j-1);
    N1 =(1/h^3)*(h^3-3*h*y3(j)^2+2*y3(j)^3);
    N2 =(1/h^2)*(y3(j)*h^2-2*h*y3(j)^2+y3(j)^3);
    N3 =(1/h^3)*(3*h*y3(j)^2-2*y3(j)^3);
    N4 =(1/h^2)*(-h*y3(j)^2+y3(j)^3);
    x3(j) =[N1 N2 N3 N4]*[0;U(2);U(6);U(8)]+L;
end
X1N2=x1;Y1N2=y1;X2N2=x2;Y2N2=y2;X3N2=x3;Y3N2=y3;
%% Second-order Linear Analysis
for ii=1:ii 
% calculate Force
 F = [P7*ii;P9;P12];
% assemble global K matrix
        K=[D+D B B;
           B A+AA C;
           B C A+AA];
        uL2 = K^-1 *F;        
% Force and Dispalcement for each increment
    F1L2(ii) = F(1);
    U1L(ii) = uL2(1);
    U2L(ii) = uL2(2);
    U3L(ii) = uL2(3);
end
D2L=[0;U1L'];F2L=[0;F1L2'];% Force-displacement of 2nd-order Linear analysis
%% First-order Nonlinear Analysis bilinear fitting
SIZE=size(D1,1);
for i=1:SIZE-1;
    hh(i) = D1(i+1)-D1(i);
    Aa(i)=(F1(i)+F1(i+1))*0.5*hh(i);
end
Area=sum(Aa);k0 =F1(50)/D1(50);
fiy = (F1(i+1)*max(D1)*0.5-Area)/(F1(i+1)*0.5 - k0*max(D1)*0.5);
Fy = k0*fiy;
X1 = [0 fiy max(D1)];Y1 = [0 Fy F1(i+1)];
%% Second-order Nonlinear Analysis bilinear fitting
SIZE=size(D2,1);
for i=1:SIZE-1;
    hhl(i) = D2(i+1)-D2(i);
    Aal(i)=(F2(i)+F2(i+1))*0.5*hhl(i);
end
Area=sum(Aal);k0 =F2(50)/D2(50);
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
fprintf('\n');
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
IMAGE=imread('PushoverAnalTwo.jpg');
image(IMAGE);axis image;axis off;
figure(4)
IMAGE=imread('PushoverAnalOfNonlinSemiRigidSpringWithElaFrameUnderDistriLoadl.jpg');
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
p1=plot([0 uc1],[0 Hc1],[0 uc2],[0 Hc2],[0 Uc1],[0 HC1],'--',[0 Uc2],[0 HC2],'--');grid on;set(p1,'LineWidth',3);
legend('1st-order anal. Hing DOF(3)','1st-order anal. Hing DOF(6)','2nd-order anal. Hing DOF(3)','2nd-order anal. Hing DOF(6)','Location','NorthEastOutside');
xlabel('Rotation (rad)');ylabel('Moment (kN.mm)');
title('Moment-Rotation diagram of Connection(Spring) During the analysis','color','b');
figure(10)
p1=plot(U8,INT1_N1_f1,'black-o',UI8,INT1_N2_f1,'blue-o');grid on;set(p1,'LineWidth',2);
legend('1st-order Nonlinear Analysis','2nd-order Nonlinear Analysis','Location','NorthEastOutside');
xlabel('Displacement (mm) [DOF (8)]');ylabel('Column Internal Axial Force (kN) [DOF (8)]');
title('Internal Axial Force-Displacement [Element (1)]  diagram of Elastic Frame with Nonlin. Semi-Rigid Connection','color','b')
figure(11)
p1=plot(D1L,F1L,'black',D1,F1,'blue',X1,Y1,'--b',D2L,F2L,'red--',D2,F2,'g',X2,Y2,'--g');grid on;set(p1,'LineWidth',3);
legend('1st-order Linear Analysis','1st-order Nonlinear Analysis','1st-order Nonlinear Curve fitted','2nd-order Linear Analysis','2nd-order Nonlinear Analysis','2nd-order Nonlinear Curve fitted','Location','NorthEastOutside');
xlabel('Displacement (mm) [DOF (7)]');ylabel('Force (kN) [DOF (7)]');
title(['Force-Disp. diagram of Elastic Frame with Nonlin. Semi-Rigid Connection',' 1-order Nonlin. Ductility Rito: ',num2str(Z3),'   2-order Nonlin. Ductility Rito: ',num2str(Z4)],'color','b')
figure(12)
p1=plot(U7,TBSN1,'black',UI7,TBSN2,'blue--');grid on;set(p1,'LineWidth',3);
legend('1st-order Nonlinear Analysis','2nd-order Nonlinear Analysis','Location','NorthEastOutside');
xlabel('Displacement (mm) [DOF (7)]');ylabel('Base Shear (kN) [DOF (1)+DOF (4)]');
title('Base Shear-Displacement diagram of Elastic Frame with Nonlin. Semi-Rigid Connection','color','b')
figure(13)
p1=plot(U7,Kz1,'black',UI7,Kz2,'blue--');grid on;set(p1,'LineWidth',3);
legend('1st-order Nonlinear Analysis','2nd-order Nonlinear Analysis','Location','NorthEastOutside');
xlabel('Displacement (mm) [DOF (7)]');ylabel('Structural Lateral Stiffness (kN/mm) [DOF (7)]');
title('Structural Lateral Stiffness-Displacement diagram of Elastic Frame with Nonlin. Semi-Rigid Connection','color','b')
figure(14)
p1=plot([U3+U6],TMSN1,'black',[UI3+UI6],TMSN2,'blue--');grid on;set(p1,'LineWidth',3);
legend('1st-order Nonlinear Analysis','2nd-order Nonlinear Analysis','Location','NorthEastOutside');
xlabel('Rotation (rad) [DOF (3)+DOF (6)]');ylabel('Base Moment (kN.mm) [DOF (3)+DOF (6)]');
title(['Base Moment-Rotation diagram of Elastic Frame with Nonlin. Semi-Rigid Connection',' 1-order Nonlin. Ductility Rito: ',num2str(Z3),'   2-order Nonlin. Ductility Rito: ',num2str(Z4)],'color','b')
figure(15)
p1=plot(corx,cory,'--b',X1N1,Y1N1,'--r',X2N1,Y2N1,'--r',X3N1,Y3N1,'--r',X1N2,Y1N2,'-.g',X2N2,Y2N2,'-.g',X3N2,Y3N2,'-.g');grid on;set(p1,'LineWidth',3);
legend('BLUE: Not Loading','RED: 1st-order Nonlinear Analysis','GREEN: 2nd-order Nonlinear Analysis','Location','NorthOutside');
xlabel('X (mm)');ylabel('Y (mm)');
title('Small Deflection Theory - Last Step Displacement of Elastic Frame with Nonlin. Semi-Rigid Connection','color','b')
figure(16)
subplot(3,1,1)
p1=plot([0 h],[INT1_N1_f1(i) -INT1_N1_f4(i)],'black',[0 h],[INT1_N2_f1(ii) -INT1_N2_f4(ii)],'b--');grid on;set(p1,'LineWidth',2);
xlabel('column length (mm)');ylabel('Axial Force (kN)');
title('Last Step Axial Force Diagram Of Element 1 - BLACK: 1st-order Nonlinear Analysis - BLUE: 2nd-order Nonlinear Analysis','color','b')
subplot(3,1,2)
p1=plot([0 L],[INT3_N1_f1(i) -INT3_N1_f4(i)],'black',[0 L],[INT3_N2_f1(ii) -INT3_N2_f4(ii)],'b--');grid on;set(p1,'LineWidth',2);
xlabel('Beam length (mm)');ylabel('Axial Force (kN)');title('Element 3','color','b');
subplot(3,1,3)
p1=plot([0 h],[INT2_N1_f1(i) -INT2_N1_f4(i)],'black',[0 h],[INT2_N2_f1(ii) -INT2_N2_f4(ii)],'b--');grid on;set(p1,'LineWidth',2);
xlabel('column length (mm)');ylabel('Axial Force (kN)');title('Element 2','color','b');
figure(17)
subplot(3,1,1)
p1=plot([0 h],[INT1_N1_f2(i) -INT1_N1_f5(i)],'black',[0 h],[INT1_N2_f2(ii) -INT1_N2_f5(ii)],'b--');grid on;set(p1,'LineWidth',2);
xlabel('column length (mm)');ylabel('Shear (kN)');
title('Last Step Shear Diagram Of Element 1 - BLACK: 1st-order Nonlinear Analysis - BLUE: 2nd-order Nonlinear Analysis','color','b')
subplot(3,1,2)
p1=plot([0 L],[INT3_N1_f2(i) -INT3_N1_f5(i)],'black',[0 L],[INT3_N2_f2(ii) -INT3_N2_f5(ii)],'b--');grid on;set(p1,'LineWidth',2);
xlabel('Beam length (mm)');ylabel('Shear (kN)');title('Element 3','color','b');
subplot(3,1,3)
p1=plot([0 h],[INT2_N1_f2(i) -INT2_N1_f5(i)],'black',[0 h],[INT2_N2_f2(ii) -INT2_N2_f5(ii)],'b--');grid on;set(p1,'LineWidth',2);
xlabel('column length (mm)');ylabel('Shear (kN)');title('Element 2','color','b');
figure(18)
subplot(3,1,1)
p1=plot([0 h],[INT1_N1_f3(1) -INT1_N1_f6(1)],'black',[0 h],[INT1_N1_f3(i) -INT1_N1_f6(i)],'--black',[0 h],[INT1_N2_f3(1) -INT1_N2_f6(1)],'b',[0 h],[INT1_N2_f3(ii) -INT1_N2_f6(ii)],'b--');grid on;set(p1,'LineWidth',2);
xlabel('column length (mm)');ylabel('Moment (kN.mm)');
title('First Step(--) Last Step(- -) Moment Diagram Of Element 1 - BLACK: 1st-order Nonlinear Analysis - BLUE: 2nd-order Nonlinear Analysis','color','b')
subplot(3,1,2)
p1=plot([0 L],[INT3_N1_f3(1) -INT3_N1_f6(1)],'black',[0 L],[INT3_N1_f3(i) -INT3_N1_f6(i)],'--black',[0 L],[INT3_N2_f3(1) -INT3_N2_f6(1)],'b',[0 L],[INT3_N2_f3(ii) -INT3_N2_f6(ii)],'b--');grid on;set(p1,'LineWidth',2);
xlabel('Beam length (mm)');ylabel('Moment (kN.mm)');title('Element 3','color','b');
subplot(3,1,3)
p1=plot([0 h],[INT2_N1_f3(1) -INT2_N1_f6(1)],'black',[0 h],[INT2_N1_f3(i) -INT2_N1_f6(i)],'--black',[0 h],[INT2_N2_f3(1) -INT2_N2_f6(1)],'b',[0 h],[INT2_N2_f3(ii) -INT2_N2_f6(ii)],'b--');grid on;set(p1,'LineWidth',2);
xlabel('column length (mm)');ylabel('Moment (kN.mm)');title('Element 2','color','b');
%%  print time of computation
totaltime = cputime - starttime;
fprintf('\nTotal time (s)= %7.4f \n\n',totaltime)
%% Output Data to .txt file
fid = fopen('PushoverAnalTwoS-OutPut.txt','w');
fprintf(fid,'======================= Result ========================\n');
fprintf(fid,'=== 1st-order Nonlinear ==+== 2nd-order Nonlinear =====\n');
fprintf(fid,'Disp.(D7) Base Shear(D1+D4) Disp.(D7) Base Shear(D1+D4)\n');
fprintf(fid,'=======================================================\n');
fprintf(fid,'  (mm)         (kN)          (mm)          (kN)\n');
fprintf(fid,'-------------------------------------------------------\n');
fprintf(fid,'%10.3f %10.3f %10.3f %10.3f\n',[X1;Y1;X2;Y2]);
fprintf(fid,'=======================================================\n');
fprintf(fid,'\n');
fprintf(fid,'+----------------------------------------------------------------------+\n');
fprintf(fid,' Semi-Rigid Column Connection Ductility Rito is (Tu/Ty): %6.3f\n',Z1);
fprintf(fid,' 1st-order Nonlinear Ductility Rito is (Du/Dy): %6.3f\n',Z3);
fprintf(fid,' 2nd-order Nonlinear Ductility Rito is (Du/Dy): %5.3f\n',Z4);
fprintf(fid,' 1st-order Nonlinear Over Strength Ratio is (Fu/Fy): %6.3f\n',Z5);
fprintf(fid,' 2nd-order Nonlinear Over Strength Ratio is (Fu/Fy): %5.3f\n',Z6);
fprintf(fid,' 1st-order Nonlinear Initial Strucural stiffness is (Ke): %6.3f [kN/mm]\n',Kyf);
fprintf(fid,' 1st-order Nonlinear Tangent Strucural stiffness is (Kt): %5.3f [kN/mm]\n',Ktf);
fprintf(fid,' 2nd-order Nonlinear Initial Strucural stiffness is (Ke): %6.3f [kN/mm]\n',Kys);
fprintf(fid,' 2nd-order Nonlinear Tangent Strucural stiffness is (Kt): %5.3f [kN/mm]\n',Kts);
fprintf(fid,'+----------------------------------------------------------------------+\n');
fprintf(fid,'\n');
fprintf(fid,'==================== Internal Moment for each degree of freedom ====================\n');
fprintf(fid,'+ ============================== 1st-order Nonlinear ============================= +\n');
fprintf(fid,'  step(i)     (f3)       (f9c)       (f9b)       (f12b)       (f6)        (f12c)  \n');
fprintf(fid,'------------------------------------------------------------------------------------\n');
fprintf(fid,'%10.0f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n',[I1;INT1_N1_f3;-INT1_N1_f6;INT3_N1_f3;-INT3_N1_f6;INT2_N1_f3;-INT2_N1_f6]);
fprintf(fid,'+ ============================== 2nd-order Nonlinear ============================= +\n');
fprintf(fid,'  step(i)     (f3)       (f9c)       (f9b)       (f12b)       (f6)        (f12c)  \n');
fprintf(fid,'------------------------------------------------------------------------------------\n');
fprintf(fid,'%10.0f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n',[I2;INT1_N2_f3;-INT1_N2_f6;INT3_N2_f3;-INT3_N2_f6;INT2_N2_f3;-INT2_N2_f6]); 
fprintf(fid,'+ ================================================================================ +\n');
fclose(fid);