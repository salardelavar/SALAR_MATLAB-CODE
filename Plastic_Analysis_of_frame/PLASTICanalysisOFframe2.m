%***********************************************************%
%              >> IN THE NAME OF GOD <<                     %
% Plastic Analysis of plane frame with Linear Programming   %
%-----------------------------------------------------------%
%     This program is written by salar delavar ghashghaei   %  
%          E-mail:salar.d.ghashghaei@gmail.com              %
%-----------------------------------------------------------%
%                     Vertical Distributed Load             %
%                                 w                         %
%                          ||||||||||||||                   %
%        Horizontal Force->|************|   -               %
%                          |     Mp     |   |               %
%                          |Mp       Mp |   h               %
%                          |            |   |               %
%                          =            =   -               %
%                          |---- L -----|                   %
%***********************************************************%
close all;clear all;clc;
disp('unit: Free')
h=6000;% Height of column
L=4000;% Length of half beam
w=10;% Vertical Distributed load(w)
H=10;% Horizontal Force
n=50; % number of Plastic moment
%Z=Plastic section Modulus --------- Mp=Plasic Moment
%% In the example Plasic Moment of beam (Mpb) and Plasic Moment of column (Mpc)and we have 1 mechansims
% by linear Programing Calculate the Horizontal Force(H) and Vertical Distributed load(w)
%%% monitor cpu time
starttime = cputime;
for i=1:n;
    disp('>>>>======================== Start ========================<<<<')
    fprintf('                    Plastic Moment-%g\n',i);
H=i*w;
z1=(3*L*w + (3*L^2*w^2 + 4*H*h*w)^(1/2))/(2*w);
z2=(3*L*w - (3*L^2*w^2 + 4*H*h*w)^(1/2))/(2*w);
z=min(z1,z2);
Mp= [(H*h+0.5*w*L*z)/(2+(L/(L-z)))]
Hwi=H/w;Hw(i)=[Hwi];Mpa(i)=[Mp];
disp('>>>>========================= End =========================<<<<')
end
Plastic_Moment____Hw=[Mpa' Hw']
%% imaging
figure (1)
IMAGE=imread('PLASTICanalysisOFframe2.jpg');
image(IMAGE);axis image;axis off;
%% Plot
figure(2)
semilogy(Hw,Mpa,'-o')
grid on;
    xlabel('Horizontal Force(H) / Vertical Distributed load(w)');
    ylabel('Plastic Moment(Mp)  (base 10 logarithm)');
    title('Interaction diagram for plastic analysis of plane frame');
%plot(HV,Mpa,'-o');
    %xlabel('Horizontal Force(H) / Vertical Force(V)');
    %ylabel('Plastic Moment(Mp)  ');
    %title('Interaction diagram for plastic analysis of plane frame')    
    %%  print time of computation
disp('--------------------');
totaltime = cputime - starttime;
fprintf('\nTotal time (s)= %7.4f \n\n',totaltime)
disp('--------------------');
