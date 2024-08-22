%***********************************************************%
%              >> IN THE NAME OF GOD <<                     %
% Plastic Analysis of plane frame with Linear Programming   %
%-----------------------------------------------------------%
%     This program is written by salar delavar ghashghaei   %  
%          E-mail:salar.d.ghashghaei@gmail.com              %
%-----------------------------------------------------------%
%                          Vertical Force                   %
%                                 |                         %
%                                 V                         %
%        Horizontal Force->|************|   -               %
%                          |     Mpb    |   |               %
%                          |Mpc      Mpc|   h               %
%                          |            |   |               %
%                          =            =   -               %
%                          |---- L -----|                   %
%***********************************************************%
close all;clear all;clc;
disp('unit: Free')
h1=6000;% Height of column1
L1=10000;% Length of half beam1
L2=10000;% Length of half beam2
H1=10;% Horizontal Force
V=100;%  Vertical Force
n=50; % number of Plastic moment
%Z=Plastic section Modulus --------- Mp=Plasic Moment
%% In the example Plasic Moment of beam (Mpb) and Plasic Moment of column (Mpc)and we have at least 6 mechansims
% by linear Programing Calculate the Horizontal Force(H) and Vertical Force(V)
%%% monitor cpu time
starttime = cputime;
for i=1:n;
disp('>>>>======================== Start ========================<<<<')
    fprintf('                    Plastic Moment-%g\n',i);
H2=i*H1;
f = [h1; L1];
%   Mpb Mpc
A = -[4 0
      0 6
      0 4
      3 1.5];    
b =-[V*0.5*L1;
     H2*4*h1;
     H1*h1+H2*h1];
lb = [1000 1000]; % lower bound
ub = [25000 25000]; % upper bound
[x,fval,exitflag,output,lambda] = linprog(f,A,b,[],[],lb,ub)
Mpb=x(1,1);Mpc=x(2,1);
HVi=H2/V;HV(i)=[HVi]';
Mpbb(i)=Mpb;Mpcc(i)=Mpc;
Mp=Mpc/Mpb;Mpa(i)=[Mp]';
disp('>>>>========================= End =========================<<<<')
end
Plastic_Moment____HV=[Mpa' HV']
%% imaging
figure (1)
IMAGE=imread('PLASTICanalysisOFframe8.jpg');
image(IMAGE);axis image;axis off;
%% Plot
figure(2)
%semilogy(HV,Mpa,'-o')
%    xlabel('Horizontal Force(H) / Vertical Force(V)');
%    ylabel('Beam Plastic Moment(Mpb) / Column Plastic Moment(Mpc)  (base 10 logarithm)');
%    title('Interaction diagram for plastic analysis of plane frame')
plot(HV,Mpa,'-o');
    xlabel('Horizontal Force(H2) / Vertical Force(V)');
    ylabel('Beam Plastic Moment(Mpb) / Column Plastic Moment(Mpc)');
    title('Interaction diagram for plastic analysis of plane frame')
grid on;
figure(3)
plot(HV,Mpbb,'-o',HV,Mpcc,'r-o');
    xlabel('Horizontal Force(H2) / Vertical Force(V)');
    ylabel('Beam Plastic Moment(Mpb) or Column Plastic Moment(Mpc)');
    title('Interaction diagram for plastic analysis of plane frame')
    legend('Beam Plastic Moment','column Plastic Moment','Location','NorthEastOutside');grid on;
disp('--------------------');
totaltime = cputime - starttime;
fprintf('\nTotal time (s)= %7.4f \n\n',totaltime)
disp('--------------------');

