%**************************************************************%
%                  >> IN THE NAME OF ALLAH <<                  %
% This program analysis Plate Element on elastic foundation    %           
% by finite difference Method                                  %
%--------------------------------------------------------------%
%      This program is written by salar delavar ghashghaei     %
%              E-mail: salar.d.ghashghaei@gmail.com            %
%--------------------------------------------------------------%
%                                                  P           %
%            -    +-------------------+            |           %
%            |    |                   |            V           %
%                 |                   |     +--------------+   %
%                 |                   |           spring       %
% Y         Ly    |         +         |     |<-   Lx     ->|   %
% ^               |                   |            P           %
% |               |                   |            |           %
% |          |    |                   |            V           %
% +-----> X  -    +-------------------+     +--------------+   %
% Z                                              spring        %
%                 |<-      Lx       ->|     |<-   Ly     ->|   %
%**************************************************************%
%             Meshing and node condition
%47       36      37       38       39       40       49
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        %       %        %        %        %        %
%29      %30     %31      %32      %33      %34      %35
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        %       %        %        %        %        %
%22      %23     %24      %25      %26      %27      %28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        %       %        %        %        %        %
%15      %16     %17      %18      %19      %20      %21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        %       %        %        %        %        %
%8       %9      %10      %11      %12      %13      %14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        %       %        %        %        %        %
%1       %2      %3       %4       %5       %6       %7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        %       %        %        %        %        %
%48     %41     %42      %43      %44      %45      %46
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;close all;clear all;
P=-100000;%input('Define P (N):');
Lx=1200;%input('Define Lx (mm):');
Ly=1200;%input('Define Ly (mm):');
h=1000;%input('Define h (mm):');
Ny=0.15;%input('Define Ny:');
E=2.1e5;%input('Define E (N/mm^2):');
k=0.15; %input('Define Soil stiffness (N/mm):')
%%% monitor cpu time
starttime = cputime;
    landaX=Lx/6;landaY=Ly/6;
    for i=1:1:7;
    X=(i-1)*landaX;Y=(i-1)*landaY;   
    x(i)=X;
    y(i)=Y;
    end
%% Finite-difference coefficients
a=(landaX/landaY)^2;N=(E*(h^3))/(12*(1-Ny));
A=6+6*a^2+8*a;B=-4*(1+a);C=-4*a*(1+a);D=2*a;E=a^2;F=1;
SoilSt=k*eye(49,49); % soil Stiffness Matrix
%% Equivalent stiffness Matrix
%1      1 2 3 4 5 6 7 8 9 10
k01_10=[A B F 0 0 0 0 C D 0];
k01_20=[0 0 0 0 E 0 0 0 0 0];
k01_30=[0 0 0 0 0 0 0 0 0 0];
k01_40=[0 0 0 0 0 0 0 0 0 0];
k01_49=[D 0 0 0 0 0 0 C 0];
k1=[k01_10 k01_20 k01_30 k01_40 k01_49];
%2      1 2 3 4 5 6 7 8 9 10
k02_10=[B A B F 0 0 0 D C D];
k02_20=[0 0 0 0 0 E 0 0 0 0];
k02_30=[0 0 0 0 0 0 0 0 0 0];
k02_40=[0 0 0 0 0 0 0 0 0 0];
k02_49=[C D 0 0 0 0 0 D 0];
k2=[k02_10 k02_20 k02_30 k02_40 k02_49];
%3      1 2 3 4 5 6 7 8 9 10
k03_10=[F B A B F 0 0 0 D C];
k03_20=[D 0 0 0 0 0 E 0 0 0];
k03_30=[0 0 0 0 0 0 0 0 0 0];
k03_40=[0 0 0 0 0 0 0 0 0 0];
k03_49=[D C D 0 0 0 0 0 0];
k3=[k03_10 k03_20 k03_30 k03_40 k03_49];
%4      1 2 3 4 5 6 7 8 9 10
k04_10=[0 F B A B F 0 0 0 D];
k04_20=[C D 0 0 0 0 0 E 0 0];
k04_30=[0 0 0 0 0 0 0 0 0 0];
k04_40=[0 0 0 0 0 0 0 0 0 0];
k04_49=[0 D C D 0 0 0 0 0];
k4=[k04_10 k04_20 k04_30 k04_40 k04_49];
%5      1 2 3 4 5 6 7 8 9 10
k05_10=[0 0 F B A B F 0 0 0];
k05_20=[D C D 0 0 0 0 0 E 0];
k05_30=[0 0 0 0 0 0 0 0 0 0];
k05_40=[0 0 0 0 0 0 0 0 0 0];
k05_49=[0 0 D C D 0 0 0 0];
k5=[k05_10 k05_20 k05_30 k05_40 k05_49];
%6      1 2 3 4 5 6 7 8 9 10
k06_10=[0 0 0 F B A B 0 0 0];
k06_20=[0 D C D 0 0 0 0 0 E];
k06_30=[0 0 0 0 0 0 0 0 0 0];
k06_40=[0 0 0 0 0 0 0 0 0 0];
k06_49=[0 0 0 D C D 0 0 0];
k6=[k06_10 k06_20 k06_30 k06_40 k06_49];
%7      1 2 3 4 5 6 7 8 9 10
k07_10=[0 0 0 0 F B A 0 0 0];
k07_20=[0 0 D C 0 0 0 0 0 0];
k07_30=[E 0 0 0 0 0 0 0 0 0];
k07_40=[0 0 0 0 0 0 0 0 0 0];
k07_49=[0 0 0 0 D C 0 0 0];
k7=[k07_10 k07_20 k07_30 k07_40 k07_49];
%8      1 2 3 4 5 6 7 8 9 10
k08_10=[C D 0 0 0 0 0 A B F];
k08_20=[0 0 0 0 C D 0 0 0 0];
k08_30=[0 E 0 0 0 0 0 0 0 0];
k08_40=[0 0 0 0 0 0 0 0 0 0];
k08_49=[0 0 0 0 0 0 0 E 0];
k8=[k08_10 k08_20 k08_30 k08_40 k08_49];
%9      1 2 3 4 5 6 7 8 9 10
k09_10=[D C D 0 0 0 0 B A B];
k09_20=[F 0 0 0 D C D 0 0 0];
k09_30=[0 0 E 0 0 0 0 0 0 0];
k09_40=[0 0 0 0 0 0 0 0 0 0];
k09_49=[E 0 0 0 0 0 0 0 0];
k9=[k09_10 k09_20 k09_30 k09_40 k09_49];
%10      1 2 3 4 5 6 7 8 9 10
k10_10= [0 D C D 0 0 0 F B A];
k10_20= [B F 0 0 0 D C D 0 0];
k10_30= [0 0 0 E 0 0 0 0 0 0];
k10_40= [0 0 0 0 0 0 0 0 0 0];
k10_49= [0 E 0 0 0 0 0 0 0];
k10=[k10_10 k10_20 k10_30 k10_40 k10_49];
%11      1 2 3 4 5 6 7 8 9 10
k11_10= [0 0 D C D 0 0 0 F B];
k11_20= [A B F 0 0 0 D C D 0];
k11_30= [0 0 0 0 E 0 0 0 0 0];
k11_40= [0 0 0 0 0 0 0 0 0 0];
k11_49= [0 0 E 0 0 0 0 0 0];
k11=[k11_10 k11_20 k11_30 k11_40 k11_49];
%12      1 2 3 4 5 6 7 8 9 10
k12_10= [0 0 0 D C D 0 0 0 F];
k12_20= [B A B F 0 0 0 D C D];
k12_30= [0 0 0 0 0 E 0 0 0 0];
k12_40= [0 0 0 0 0 0 0 0 0 0];
k12_49= [0 0 0 E 0 0 0 0 0];
k12=[k12_10 k12_20 k12_30 k12_40 k12_49];
%13      1 2 3 4 5 6 7 8 9 10
k13_10= [0 0 0 0 D C D 0 0 0];
k13_20= [F B A B 0 0 0 0 D C];
k13_30= [D 0 0 0 0 0 E 0 0 0];
k13_40= [0 0 0 0 0 0 0 0 0 0];
k13_49= [0 0 0 0 E 0 0 0 0];
k13=[k13_10 k13_20 k13_30 k13_40 k13_49];
%14      1 2 3 4 5 6 7 8 9 10
k14_10= [0 0 0 0 0 D C 0 0 0];
k14_20= [0 F B A 0 0 0 0 0 D];
k14_30= [C 0 0 0 0 0 0 E 0 0];
k14_40= [0 0 0 0 0 0 0 0 0 0];
k14_49= [0 0 0 0 0 E 0 0 0];
k14=[k14_10 k14_20 k14_30 k14_40 k14_49];
%15      1 2 3 4 5 6 7 8 9 10
k15_10= [E 0 0 0 0 0 0 C D 0];
k15_20= [0 0 0 0 A B F 0 0 0];
k15_30= [0 C D 0 0 0 0 0 E 0];
k15_40= [0 0 0 0 0 0 0 0 0 0];
k15_49= [0 0 0 0 0 0 0 0 0];
k15=[k15_10 k15_20 k15_30 k15_40 k15_49];
%16      1 2 3 4 5 6 7 8 9 10
k16_10= [0 E 0 0 0 0 0 D C D];
k16_20= [0 0 0 0 B A B F 0 0];
k16_30= [0 D C D 0 0 0 0 0 E];
k16_40= [0 0 0 0 0 0 0 0 0 0];
k16_49= [0 0 0 0 0 0 0 0 0];
k16=[k16_10 k16_20 k16_30 k16_40 k16_49];
%17      1 2 3 4 5 6 7 8 9 10
k17_10= [0 0 E 0 0 0 0 0 D C];
k17_20= [D 0 0 0 F B A B F 0];
k17_30= [0 0 D C D 0 0 0 0 0];
k17_40= [E 0 0 0 0 0 0 0 0 0];
k17_49= [0 0 0 0 0 0 0 0 0];
k17=[k17_10 k17_20 k17_30 k17_40 k17_49];
%18      1 2 3 4 5 6 7 8 9 10
k18_10= [0 0 0 E 0 0 0 0 0 D];
k18_20= [C D 0 0 0 F B A B F];
k18_30= [0 0 0 D C D 0 0 0 0];
k18_40= [0 E 0 0 0 0 0 0 0 0];
k18_49= [0 0 0 0 0 0 0 0 0];
k18=[k18_10 k18_20 k18_30 k18_40 k18_49];
%19      1 2 3 4 5 6 7 8 9 10
k19_10= [0 0 0 0 E 0 0 0 0 0];
k19_20= [D C D 0 0 0 F B A B];
k19_30= [F 0 0 0 D C D 0 0 0];
k19_40= [0 0 E 0 0 0 0 0 0 0];
k19_49= [0 0 0 0 0 0 0 0 0];
k19=[k19_10 k19_20 k19_30 k19_40 k19_49];
%20      1 2 3 4 5 6 7 8 9 10
k20_10= [0 0 0 0 0 E 0 0 0 0];
k20_20= [0 D C D 0 0 0 F B A];
k20_30= [B 0 0 0 0 D C D 0 0];
k20_40= [0 0 0 E 0 0 0 0 0 0];
k20_49= [0 0 0 0 0 0 0 0 0];
k20=[k20_10 k20_20 k20_30 k20_40 k20_49];
%21      1 2 3 4 5 6 7 8 9 10
k21_10= [0 0 0 0 0 0 E 0 0 0];
k21_20= [0 0 D C 0 0 0 0 F B];
k21_30= [A 0 0 0 0 0 D C 0 0];
k21_40= [0 0 0 0 E 0 0 0 0 0];
k21_49= [0 0 0 0 0 0 0 0 0];
k21=[k21_10 k21_20 k21_30 k21_40 k21_49];
%22      1 2 3 4 5 6 7 8 9 10
k22_10= [0 0 0 0 0 0 0 E 0 0];
k22_20= [0 0 0 0 C D 0 0 0 0];
k22_30= [0 A B F 0 0 0 0 C D];
k22_40= [0 0 0 0 0 0 0 0 0 0];
k22_49= [0 0 0 0 0 0 E 0 0];
k22=[k22_10 k22_20 k22_30 k22_40 k22_49];
%23      1 2 3 4 5 6 7 8 9 10
k23_10= [0 0 0 0 0 0 0 0 E 0];
k23_20= [0 0 0 0 D C D 0 0 0];
k23_30= [0 B A B F 0 0 0 D C];
k23_40= [D 0 0 0 0 E 0 0 0 0];
k23_49= [0 0 0 0 0 0 0 0 0];
k23=[k23_10 k23_20 k23_30 k23_40 k23_49];
%24      1 2 3 4 5 6 7 8 9 10
k24_10= [0 0 0 0 0 0 0 0 0 E];
k24_20= [0 0 0 0 0 D C D 0 0];
k24_30= [0 F B A B F 0 0 0 D];
k24_40= [C D 0 0 0 0 E 0 0 0];
k24_49= [0 0 0 0 0 0 0 0 0];
k24=[k24_10 k24_20 k24_30 k24_40 k24_49];
%25      1 2 3 4 5 6 7 8 9 10
k25_10= [0 0 0 0 0 0 0 0 0 0];
k25_20= [E 0 0 0 0 0 D C D 0];
k25_30= [0 0 F B A B F 0 0 0];
k25_40= [D C D 0 0 0 0 E 0 0];
k25_49= [0 0 0 0 0 0 0 0 0];
k25=[k25_10 k25_20 k25_30 k25_40 k25_49];
%26      1 2 3 4 5 6 7 8 9 10
k26_10= [0 0 0 0 0 0 0 0 0 0];
k26_20= [0 E 0 0 0 0 0 D C D];
k26_30= [0 0 0 F B A B F 0 0];
k26_40= [0 D C D 0 0 0 0 E 0];
k26_49= [0 0 0 0 0 0 0 0 0];
k26=[k26_10 k26_20 k26_30 k26_40 k26_49];
%27      1 2 3 4 5 6 7 8 9 10
k27_10= [0 0 0 0 0 0 0 0 0 0];
k27_20= [0 0 E 0 0 0 0 0 D C];
k27_30= [D 0 0 0 F B A B 0 0];
k27_40= [0 0 D C D 0 0 0 0 E];
k27_49= [0 0 0 0 0 0 0 0 0];
k27=[k27_10 k27_20 k27_30 k27_40 k27_49];
%28      1 2 3 4 5 6 7 8 9 10
k28_10= [0 0 0 0 0 0 0 0 0 0];
k28_20= [0 0 0 E 0 0 0 0 0 D];
k28_30= [C 0 0 0 0 F B A 0 0];
k28_40= [0 0 0 D C 0 0 0 0 0];
k28_49= [0 0 0 0 0 0 0 0 E];
k28=[k28_10 k28_20 k28_30 k28_40 k28_49];
%29      1 2 3 4 5 6 7 8 9 10
k29_10= [0 0 0 0 0 0 0 0 0 0];
k29_20= [0 0 0 0 E 0 0 0 0 0];
k29_30= [0 C D 0 0 0 0 0 A B];
k29_40= [F 0 0 0 0 D 0 0 0 0];
k29_49= [0 0 0 0 0 0 C 0 0];
k29=[k29_10 k29_20 k29_30 k29_40 k29_49];
%30      1 2 3 4 5 6 7 8 9 10
k30_10= [0 0 0 0 0 0 0 0 0 0];
k30_20= [0 0 0 0 0 E 0 0 0 0];
k30_30= [0 D C D 0 0 0 0 B A];
k30_40= [B F 0 0 0 C D 0 0 0];
k30_49= [0 0 0 0 0 0 D 0 0];
k30=[k30_10 k30_20 k30_30 k30_40 k30_49];
%31      1 2 3 4 5 6 7 8 9 10
k31_10= [0 0 0 0 0 0 0 0 0 0];
k31_20= [0 0 0 0 0 0 E 0 0 0];
k31_30= [0 0 D C D 0 0 0 F B];
k31_40= [A B F 0 0 D C D 0 0];
k31_49= [0 0 0 0 0 0 0 0 0];
k31=[k31_10 k31_20 k31_30 k31_40 k31_49];
%32      1 2 3 4 5 6 7 8 9 10
k32_10= [0 0 0 0 0 0 0 0 0 0];
k32_20= [0 0 0 0 0 0 0 E 0 0];
k32_30= [0 0 0 D C D 0 0 0 F];
k32_40= [B A B F 0 0 D C D 0];
k32_49= [0 0 0 0 0 0 0 0 0];
k32=[k32_10 k32_20 k32_30 k32_40 k32_49];
%33      1 2 3 4 5 6 7 8 9 10
k33_10= [0 0 0 0 0 0 0 0 0 0];
k33_20= [0 0 0 0 0 0 0 0 E 0];
k33_30= [0 0 0 0 D C D 0 0 0];
k33_40= [F B A B F 0 0 D C D];
k33_49= [0 0 0 0 0 0 0 0 0];
k33=[k33_10 k33_20 k33_30 k33_40 k33_49];
%34      1 2 3 4 5 6 7 8 9 10
k34_10= [0 0 0 0 0 0 0 0 0 0];
k34_20= [0 0 0 0 0 0 0 0 0 E];
k34_30= [0 0 0 0 0 D C D 0 0];
k34_40= [0 F B A B 0 0 0 D C];
k34_49= [0 0 0 0 0 0 0 0 D];
k34=[k34_10 k34_20 k34_30 k34_40 k34_49];
%35      1 2 3 4 5 6 7 8 9 10
k35_10= [0 0 0 0 0 0 0 0 0 0];
k35_20= [0 0 0 0 0 0 0 0 0 0];
k35_30= [E 0 0 0 0 0 D C 0 0];
k35_40= [0 0 F B A 0 0 0 0 D];
k35_49= [0 0 0 0 0 0 0 0 C];
k35=[k35_10 k35_20 k35_30 k35_40 k35_49];
%36      1 2 3 4 5 6 7 8 9 10
k36_10= [0 0 0 0 0 0 0 0 0 0];
k36_20= [0 0 0 0 0 0 0 0 0 0];
k36_30= [0 0 E 0 0 0 0 0 D C];
k36_40= [D 0 0 0 0 A B F 0 0];
k36_49= [0 0 0 0 0 0 B 0 0];
k36=[k36_10 k36_20 k36_30 k36_40 k36_49];
%37      1 2 3 4 5 6 7 8 9 10
k37_10= [0 0 0 0 0 0 0 0 0 0];
k37_20= [0 0 0 0 0 0 0 0 0 0];
k37_30= [0 0 0 E 0 0 0 0 0 D];
k37_40= [C D 0 0 0 B A B F 0];
k37_49= [0 0 0 0 0 0 F 0 0];
k37=[k37_10 k37_20 k37_30 k37_40 k37_49];
%38      1 2 3 4 5 6 7 8 9 10
k38_10= [0 0 0 0 0 0 0 0 0 0];
k38_20= [0 0 0 0 0 0 0 0 0 0];
k38_30= [0 0 0 0 E 0 0 0 0 0];
k38_40= [D C D 0 0 F B A B F];
k38_49= [0 0 0 0 0 0 0 0 0];
k38=[k38_10 k38_20 k38_30 k38_40 k38_49];
%39      1 2 3 4 5 6 7 8 9 10
k39_10= [0 0 0 0 0 0 0 0 0 0];
k39_20= [0 0 0 0 0 0 0 0 0 0];
k39_30= [0 0 0 0 0 E 0 0 0 0];
k39_40= [0 D C D 0 0 F B A B];
k39_49= [0 0 0 0 0 0 0 0 F];
k39=[k39_10 k39_20 k39_30 k39_40 k39_49];
%40      1 2 3 4 5 6 7 8 9 10
k40_10= [0 0 0 0 0 0 0 0 0 0];
k40_20= [0 0 0 0 0 0 0 0 0 0];
k40_30= [0 0 0 0 0 0 E 0 0 0];
k40_40= [0 0 D C D 0 0 F B A];
k40_49= [0 0 0 0 0 0 0 0 0];
k40=[k40_10 k40_20 k40_30 k40_40 k40_49];
%41      1 2 3 4 5 6 7 8 9 10
k41_10= [D C D 0 0 0 0 0 E 0];
k41_20= [0 0 0 0 0 0 0 0 0 0];
k41_30= [0 0 0 0 0 0 0 0 0 0];
k41_40= [0 0 0 0 0 0 0 0 0 0];
k41_49= [A B F 0 0 0 0 0 0];
k41=[k41_10 k41_20 k41_30 k41_40 k41_49];
%42      1 2 3 4 5 6 7 8 9 10
k42_10= [0 D C D 0 0 0 0 0 E];
k42_20= [0 0 0 0 0 0 0 0 0 0];
k42_30= [0 0 0 0 0 0 0 0 0 0];
k42_40= [0 0 0 0 0 0 0 0 0 0];
k42_49= [B A B F 0 0 0 F 0];
k42=[k42_10 k42_20 k42_30 k42_40 k42_49];
%43      1 2 3 4 5 6 7 8 9 10
k43_10= [0 0 D C D 0 0 0 0 0];
k43_20= [E 0 0 0 0 0 0 0 0 0];
k43_30= [0 0 0 0 0 0 0 0 0 0];
k43_40= [0 0 0 0 0 0 0 0 0 0];
k43_49= [F B A B F 0 0 0 0];
k43=[k43_10 k43_20 k43_30 k43_40 k43_49];
%44      1 2 3 4 5 6 7 8 9 10
k44_10= [0 0 0 D C D 0 0 0 0];
k44_20= [0 E 0 0 0 0 0 0 0 0];
k44_30= [0 0 0 0 0 0 0 0 0 0];
k44_40= [0 0 0 0 0 0 0 0 0 0];
k44_49= [0 F B A B F 0 0 0];
k44=[k44_10 k44_20 k44_30 k44_40 k44_49];
%45      1 2 3 4 5 6 7 8 9 10
k45_10= [0 0 0 0 D C D 0 0 0];
k45_20= [0 0 E 0 0 0 0 0 0 0];
k45_30= [0 0 0 0 0 0 0 0 0 0];
k45_40= [0 0 0 0 0 0 0 0 0 0];
k45_49= [0 0 F B A B 0 0 0];
k45=[k45_10 k45_20 k45_30 k45_40 k45_49];
%46      1 2 3 4 5 6 7 8 9 10
k46_10= [0 0 0 0 0 D C 0 0 0];
k46_20= [0 0 0 E 0 0 0 0 0 0];
k46_30= [0 0 0 0 0 0 0 0 0 0];
k46_40= [0 0 0 0 0 0 0 0 0 0];
k46_49= [0 0 0 F B A 0 0 0];
k46=[k46_10 k46_20 k46_30 k46_40 k46_49];
%47      1 2 3 4 5 6 7 8 9 10
k47_10= [0 0 0 0 0 0 0 0 0 0];
k47_20= [0 0 0 0 0 0 0 0 0 0];
k47_30= [0 E 0 0 0 0 0 0 C D];
k47_40= [0 0 0 0 0 B F 0 0 0];
k47_49= [0 0 0 0 0 0 A 0 0];
k47=[k47_10 k47_20 k47_30 k47_40 k47_49];
%48      1 2 3 4 5 6 7 8 9 10
k48_10= [C D 0 0 0 0 0 E 0 0];
k48_20= [0 0 0 0 0 0 0 0 0 0];
k48_30= [0 0 0 0 0 0 0 0 0 0];
k48_40= [0 0 0 0 0 0 0 0 0 0];
k48_49= [B F 0 0 0 0 0 A 0];
k48=[k48_10 k48_20 k48_30 k48_40 k48_49];
%49      1 2 3 4 5 6 7 8 9 10
k49_10= [0 0 0 0 0 0 0 0 0 0];
k49_20= [0 0 0 0 0 0 0 0 0 0];
k49_30= [0 0 0 0 0 0 0 E 0 0];
k49_40= [0 0 0 D C 0 0 0 F B];
k49_49= [0 0 0 0 0 0 0 0 A];
k49=[k49_10 k49_20 k49_30 k48_40 k49_49];
a=(N/(landaX*landaY));
K=a*[k1;k2;k3;k4;k5;k6;k7;k8;k9;k10;k11;k12;k13;k14;k15;k16;k17;k18;
    k19;k20;k21;k22;k23;k24;k25;k26;k27;k28;k29;k30;k31;k32;k33;k34;
    k35;k36;k37;k38;k39;k40;k41;k42;k43;k44;k45;k46;k47;k48;k49]+SoilSt;
F=(P*landaX*landaY)*[0;0;0;0;0;0;0;
                     0;0;0;0;0;0;0;
                     0;0;0;1;0;0;0;
                     0;0;0;0;0;0;0;
                     0;0;0;0;0;0;0;
                     0;0;0;0;0;0;0;
                     0;0;0;0;0;0;0];
if det(K)<0
    disp('-------------------------------------------------------------------')
    disp(' Stiffness matrix is not positive definite. Zero value encountered')
    disp('-------------------------------------------------------------------')
    break
end
%% Displacement
D=inv(K)*F;
o1x=0;o2x=x(2);o3x=x(3);o4x=x(4);o5x=x(5);o6x=x(6);o7x=x(7);
o8x=0;o9x=x(2);o10x=x(3);o11x=x(4);o12x=x(5);o13x=x(6);o14x=x(7);
o15x=0;o16x=x(2);o17x=x(3);o18x=x(4);o19x=x(5);o20x=x(6);o21x=x(7);
o22x=0;o23x=x(2);o24x=x(3);o25x=x(4);o26x=x(5);o27x=x(6);o28x=x(7);
o29x=0;o30x=x(2);o31x=x(3);o32x=x(4);o33x=x(5);o34x=x(6);o35x=x(7);
o36x=0;o37x=x(2);o38x=x(3);o39x=x(4);o40x=x(5);o41x=x(6);o42x=x(7);
o43x=0;o44x=x(2);o45x=x(3);o46x=x(4);o47x=x(5);o48x=x(6);o49x=x(7);
X=[o1x;o2x;o3x;o4x;o5x;o6x;o7x;o8x;o9x;o10x;o11x;o12x;o13x;o14x;o15x;o16x;
    o17x;o18x;o19x;o20x;o21x;o22x;o23x;o24x;o25x;o26x;o27x;o28x;o29x;o30x;o31x;o32x;
    o33x;o34x;o35x;o36x;o37x;o38x;o39x;o40x;o41x;o42x;o43x;o44x;o45x;o46x;o47x;o48x;o49x];
o1y=0;o2y=y(1);o3y=y(1);o4y=y(1);o5y=y(1);o6y=y(1);o7y=y(1);
o8y=y(2);o9y=y(2);o10y=y(2);o11y=y(2);o12y=y(2);o13y=y(2);o14y=y(2);
o15y=y(3);o16y=y(3);o17y=y(3);o18y=y(3);o19y=y(3);o20y=y(3);o21y=y(3);
o22y=y(4);o23y=y(4);o24y=y(4);o25y=y(4);o26y=y(4);o27y=y(4);o28y=y(4);
o29y=y(5);o30y=y(5);o31y=y(5);o32y=y(5);o33y=y(5);o34y=y(5);o35y=y(5);
o36y=y(6);o37y=y(6);o38y=y(6);o39y=y(6);o40y=y(6);o41y=y(6);o42y=y(6);
o43y=y(7);o44y=y(7);o45y=y(7);o46y=y(7);o47y=y(7);o48y=y(7);o49y=y(7);
Y=[o1y;o2y;o3y;o4y;o5y;o6y;o7y;o8y;o9y;o10y;o11y;o12y;o13y;o14y;o15y;o16y;
    o17y;o18y;o19y;o20y;o21y;o22y;o23y;o24y;o25y;o26y;o27y;o28y;o29y;o30y;o31y;o32y;
    o33y;o34y;o35y;o36y;o37y;o38y;o39y;o40y;o41y;o42y;o43y;o44y;o45y;o46y;o47y;o48y;o49y];
Z=[D(48);D(41);D(42);D(43);D(44);D(45);D(46);D(1);D(2);D(3);D(4);D(5);D(6);D(7);D(8);D(9);D(10);D(11);D(12);D(13);D(14);
D(15);D(16);D(17);D(18);D(19);D(20);D(21);D(22);D(23);D(24);D(25);D(26);D(27);D(28);
D(29);D(30);D(31);D(32);D(33);D(34);D(35);D(47);D(36);D(37);D(38);D(39);D(40);D(49)];
   
X1=[X(1);X(2);X(3);X(4);X(5);X(6);X(7)];Y1=[Y(1);Y(2);Y(3);Y(4);Y(5);Y(6);Y(7)];Z1=[Z(1);Z(2);Z(3);Z(4);Z(5);Z(6);Z(7)];
X2=[X(8);X(9);X(10);X(11);X(12);X(13);X(14)];Y2=[Y(8);Y(9);Y(10);Y(11);Y(12);Y(13);Y(14)];Z2=[Z(8);Z(9);Z(10);Z(11);Z(12);Z(13);Z(14)];
X3=[X(15);X(16);X(17);X(18);X(19);X(20);X(21)];Y3=[Y(15);Y(16);Y(17);Y(18);Y(19);Y(20);Y(21)];Z3=[Z(15);Z(16);Z(17);Z(18);Z(19);Z(20);Z(21)];
X4=[X(22);X(23);X(24);X(25);X(26);X(27);X(28)];Y4=[Y(22);Y(23);Y(24);Y(25);Y(26);Y(27);Y(28)];Z4=[Z(22);Z(23);Z(24);Z(25);Z(26);Z(27);Z(28)];
X5=[X(29);X(30);X(31);X(32);X(33);X(34);X(35)];Y5=[Y(29);Y(30);Y(31);Y(32);Y(33);Y(34);Y(35)];Z5=[Z(29);Z(30);Z(31);Z(32);Z(33);Z(34);Z(35)];
X6=[X(36);X(37);X(38);X(39);X(40);X(41);X(42)];Y6=[Y(36);Y(37);Y(38);Y(39);Y(40);Y(41);Y(42)];Z6=[Z(36);Z(37);Z(38);Z(39);Z(40);Z(41);Z(42)];
X7=[X(43);X(44);X(45);X(46);X(47);X(48);X(49)];Y7=[Y(43);Y(44);Y(45);Y(46);Y(47);Y(48);Y(49)];Z7=[Z(43);Z(44);Z(45);Z(46);Z(47);Z(48);Z(49)];
X8=[X(1);X(8);X(15);X(22);X(29);X(36);X(43)];Y8=[Y(1);Y(8);Y(15);Y(22);Y(29);Y(36);Y(43)];Z8=[Z(1);Z(8);Z(15);Z(22);Z(29);Z(36);Z(43)];
X9=[X(2);X(9);X(16);X(23);X(30);X(37);X(44)];Y9=[Y(2);Y(9);Y(16);Y(23);Y(30);Y(37);Y(44)];Z9=[Z(2);Z(9);Z(16);Z(23);Z(30);Z(37);Z(44)];
X10=[X(3);X(10);X(17);X(24);X(31);X(38);X(45)];Y10=[Y(3);Y(10);Y(17);Y(24);Y(31);Y(38);Y(45)];Z10=[Z(3);Z(10);Z(17);Z(24);Z(31);Z(38);Z(45)];
X11=[X(4);X(11);X(18);X(25);X(32);X(39);X(46)];Y11=[Y(4);Y(11);Y(18);Y(25);Y(32);Y(39);Y(46)];Z11=[Z(4);Z(11);Z(18);Z(25);Z(32);Z(39);Z(46)];
X12=[X(5);X(12);X(19);X(26);X(33);X(40);X(47)];Y12=[Y(5);Y(12);Y(19);Y(26);Y(33);Y(40);Y(47)];Z12=[Z(5);Z(12);Z(19);Z(26);Z(33);Z(40);Z(47)];
X13=[X(6);X(13);X(20);X(27);X(34);X(41);X(48)];Y13=[Y(6);Y(13);Y(20);Y(27);Y(34);Y(41);Y(48)];Z13=[Z(6);Z(13);Z(20);Z(27);Z(34);Z(41);Z(48)];
X14=[X(7);X(14);X(21);X(28);X(35);X(42);X(49)];Y14=[Y(7);Y(14);Y(21);Y(28);Y(35);Y(42);Y(49)];Z14=[Z(7);Z(14);Z(21);Z(28);Z(35);Z(42);Z(49)];
%% imaging
figure (1)
IMAGE=imread('AnalysisPlateOnElasticFoundation.jpg');
image(IMAGE);axis image;axis off;
%% Ploting
figure (2)
P1=plot3(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,X5,Y5,Z5,X6,Y6,Z6,X7,Y7,Z7,X8,Y8,Z8,X9,Y9,Z9,X10,Y10,Z10,X11,Y11,Z11,X12,Y12,Z12,X13,Y13,Z13,X14,Y14,Z14);
xlabel('Lx (mm)');ylabel('Ly (mm)');zlabel('Displacement (mm)');set(P1,'LineWidth',5);grid on;
title('Deflection diagram of Plate on elastic foundation by finite difference method','color','b');
%% -------------------------------------------------------------------
disp('************ Result *************');
disp('===============================');
disp('       X       Y      Delta(w)');
disp('===============================');
disp([X Y Z]);
disp('===============================');
%%%  print time of computation
totaltime = cputime - starttime;
fprintf('\nMax Displacement (mm)= %7.4f \n\n',roundn(min(D),-2))
fprintf('\nTotal time (s)= %7.4f \n\n',totaltime)
disp('*********************************');
