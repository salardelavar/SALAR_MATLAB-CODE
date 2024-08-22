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
%                          |--- 2L -----|                   %
%***********************************************************%
close all;clear all;clc;
disp('unit: Free')
h=6000;% Height of column
L=4000;% Length of half beam
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
H=i*V;
f = [2*h; L];
%   Mpb Mpc
A = -[0 4
      2 2
      2 2
      4 0
      4 2
      2 4];    
b =-[V*L;
     V*L;
     H*h;
     H*h;
     H*h+V*L;
     H*h+V*L];
lb = zeros(2,1);
[x,fval,exitflag,output,lambda] = linprog(f,A,b,[],[],lb)
Mpb=x(1,1);Mpc=x(2,1);
HVi=H/V;HV(i)=[HVi]';
Mpbb(i)=Mpb;Mpcc(i)=Mpc;
Mp=Mpc/Mpb;Mpa(i)=[Mp]';
disp('>>>>========================= End =========================<<<<')
end
Plastic_Moment____HV=[Mpa' HV']
%% imaging
figure (1)
IMAGE=imread('PLASTICanalysisOFframe1.jpg');
image(IMAGE);axis image;axis off;
%% Plot
figure(2)
%semilogy(HV,Mpa,'-o')
%    xlabel('Horizontal Force(H) / Vertical Force(V)');
%    ylabel('Beam Plastic Moment(Mpb) / Column Plastic Moment(Mpc)  (base 10 logarithm)');
%    title('Interaction diagram for plastic analysis of plane frame')
plot(HV,Mpa,'-o');
    xlabel('Horizontal Force(H) / Vertical Force(V)');
    ylabel('Beam Plastic Moment(Mpb) / Column Plastic Moment(Mpc)');
    title('Interaction diagram for plastic analysis of plane frame')
grid on;
figure(3)
plot(HV,Mpbb,'-o',HV,Mpcc,'r-o');
    xlabel('Horizontal Force(H) / Vertical Force(V)');
    ylabel('Beam Plastic Moment(Mpb) or Column Plastic Moment(Mpc)');
    title('Interaction diagram for plastic analysis of plane frame')
    legend('Beam Plastic Moment','column Plastic Moment','Location','NorthEastOutside');grid on;
disp('--------------------');
totaltime = cputime - starttime;
fprintf('\nTotal time (s)= %7.4f \n\n',totaltime)
disp('--------------------');

