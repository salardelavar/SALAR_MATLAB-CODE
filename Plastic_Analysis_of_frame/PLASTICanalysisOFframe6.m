%***********************************************************%
%                >> IN THE NAME OF GOD <<                   %
% Plastic Analysis of plane frame with Linear Programming   %
%-----------------------------------------------------------%
%     This program is written by salar delavar ghashghaei   %  
%          E-mail:salar.d.ghashghaei@gmail.com              %
%-----------------------------------------------------------%
%                                 V3          V4            %
%                                 |            |            %
%                     2H ->|************|************|   -  %
%                          |     Mpb    |     Mpb    |   |  %
%                          |Mpc   V1 Mpc|Mpc   V2 Mpc|   h2 %
%                          |      |     |      |     |   |  %
%                      H ->|************|************|   -  %
%                          |     Mpb    |     Mpb    |   |  %
%                          |Mpc      Mpc|Mpc      Mpc|   h1 %
%                          |            |            |   |  %
%                          =            =            =   -  %
%                          |--- L1 -----|--- L2 -----|      %
%***********************************************************%
close all;clear all;clc;
disp('unit: Free')
h1=6000;% Height of column1
h2=6000;% Height of column2
L1=4000;% Length of half beam1
L2=4000;% Length of half beam1
V1=100;%  Vertical Force
n=100; % number of Plastic moment
%Z=Plastic section Modulus --------- Mp=Plasic Moment
%% In the example Plasic Moment of beam (Mpb) and Plasic Moment of column (Mpc)and we have at least 6 mechansims
% by linear Programing Calculate the Horizontal Force(H) and Vertical Force(V)
%%% monitor cpu time
starttime = cputime;
for i=1:n;
disp('>>>>======================== Start ========================<<<<')
    fprintf('                    Plastic Moment-%g\n',i);
H1=i*V1;H2=2*H1;
V2=2*V1;V3=3*V1;V4=4*V1;
f = [4*h1; 2*L1];
%   Mpb Mpc
A = -[4 0
      4 0
      4 0
      4 0
      0 6
      0 6
     12 6
     12 4];    
b =-[V1*0.5*L1;
     V2*0.5*L2;
     V3*0.5*L1;
     V4*0.5*L2;
     H2*h2;
     (H1+H2)*h2;
     H1*h1+H2*(h1+h2)+V1*.5*L1+V2*.5*L2+V3*.5*L1+V4*.5*L2;
     H2*h2+V1*.5*L1+V2*.5*L2+V3*.5*L1+V4*.5*L2;];
lb = zeros(2,1);
[x,fval,exitflag,output,lambda] = linprog(f,A,b,[],[],lb)
Mpb=x(1,1);Mpc=x(2,1);
HVi=H1/V1;HV(i)=[HVi]';
Mpbb(i)=Mpb;Mpcc(i)=Mpc;
Mp=Mpc/Mpb;Mpa(i)=[Mp]';
disp('>>>>========================= End =========================<<<<')
end
Plastic_Moment____HV=[Mpa' HV']
%% imaging
figure (1)
IMAGE=imread('PLASTICanalysisOFframe6.jpg');
image(IMAGE);axis image;axis off;
%% Plot
figure(2)
%semilogy(HV,Mpa,'-o')
    %xlabel('Horizontal Force(H) / Vertical Force(V)');grid on;
    %ylabel('Beam Plastic Moment(Mpb) / Column Plastic Moment(Mpc)  (base 10 logarithm)');
    %title('Interaction diagram for plastic analysis of plane frame');
plot(HV,Mpa,'-o');grid on;
    xlabel('Horizontal Force(H) / Vertical Force(V)');
    ylabel('Beam Plastic Moment(Mpb) / Column Plastic Moment(Mpc)');
    title('Interaction diagram for plastic analysis of plane frame');
figure(3)
plot(HV,Mpbb,'-o',HV,Mpcc,'r-o');
    xlabel('Horizontal Force(H) / Vertical Force(V)');
    ylabel('Beam Plastic Moment(Mpb) or Column Plastic Moment(Mpc)');
    title('Interaction diagram for plastic analysis of plane frame')
    legend('Beam Plastic Moment','column Plastic Moment','Location','NorthEastOutside');grid on;     
    %%  print time of computation
disp('--------------------');
totaltime = cputime - starttime;
fprintf('\nTotal time (s)= %7.4f \n\n',totaltime)
disp('--------------------');

