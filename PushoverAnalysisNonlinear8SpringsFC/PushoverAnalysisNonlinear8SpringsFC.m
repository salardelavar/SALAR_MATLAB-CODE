%***********************************************************%
%                  >> IN THE NAME OF GOD <<                 %
% Pushover Analysis of Nonlinear Springs with Force Control %
%-----------------------------------------------------------%
%     This program is written by salar delavar ghashghaei   %  
%            E-mail:salar.d.ghashghaei@gmail.com            %
%             Publication Date : 25 - May - 2017            %
%***********************************************************%
clear all;close all;clc
% Define Parameters in Free Unit
P = 1;
m = 10000; % number of calculation
itermax = 500;% maximum number of iterations
tolerance = 1e-12; % specified tolerance for convergence
u = 0;% initial guess value
%%% monitor cpu time
starttime = cputime;
%% Spring Properties
%      D1 F1 D2 F2 D3 F3 D4 F4 
DATA = [5 20 35 30 50 32 80 35;
        6 20 35 30 50 32 80 45;
        8 20 35 30 50 32 80 50;
        10 20 35 30 50 32 80 55;
        12 20 35 30 50 32 80 60;
        14 20 35 30 50 32 80 65;
        16 20 35 30 50 32 80 35;
        9 20 35 30 50 32 80 35];
        
Dmax = max(DATA(:,5)); % [mm] Max displacement 
for i=1:8
Rk1(i)=(DATA(i,2)-0)/(DATA(i,1)-0);
Rk2(i)=(DATA(i,4)-DATA(i,2))/(DATA(i,3)-DATA(i,1));
Rk3(i)=(DATA(i,6)-DATA(i,4))/(DATA(i,5)-DATA(i,3));
Rk4(i)=(DATA(i,8)-DATA(i,6))/(DATA(i,7)-DATA(i,5));
end
f=zeros(8,1);
%% Nonlinear Springs Analysis
disp('#################################################');
disp('#    Pushover Analysis of Nonlinear Springs     #');
disp('#################################################');
% Gradually increase the applied load
for i=1:m
    % Define the applied load
        F = [P*i];
        Kini=0;
        for (j=1:8)
        if and(abs(f(j)) >= 0,abs(f(j)) <= DATA(j,2))
            K(j) = Rk1(j);
        elseif and(abs(f(j)) > DATA(j,2),abs(f(j))<= DATA(j,4))
            K(j) = (DATA(j,2)+Rk2(j)*(abs(u)-DATA(j,1)))/abs(u);
        elseif and(abs(f(j)) > DATA(j,4),abs(f(j))<= DATA(j,6))
            K(j) = (DATA(j,4)+Rk3(j)*(abs(u)-DATA(j,3)))/abs(u);
        elseif and(abs(f(j)) > DATA(j,6),abs(f(j)) <= DATA(j,8))
            K(j) = (DATA(j,6)+Rk4(j)*(abs(u)-DATA(j,5)))/abs(u);
        else 
            K(j) = 0;    
            
        end
        Kini = Kini+K(j);
        end
        
        it = 0; % initialize iteration count
        residual = 100; % initialize residual
        while (residual > tolerance)
        Ko=0;
        for (j=1:8)
        if and(abs(f(j)) >= 0,abs(f(j)) <= DATA(j,2))
            Kt(j) = Rk1(j);
        elseif and(abs(f(j)) > DATA(j,2),abs(f(j))<= DATA(j,4))
            Kt(j) = (DATA(j,2)+Rk2(j)*(abs(u)-DATA(j,1)))/abs(u);
        elseif and(abs(f(j)) > DATA(j,4),abs(f(j))<= DATA(j,6))
            Kt(j) = (DATA(j,4)+Rk3(j)*(abs(u)-DATA(j,3)))/abs(u);
        elseif and(abs(f(j)) > DATA(j,6),abs(f(j)) <= DATA(j,8))
            Kt(j) = (DATA(j,6)+Rk4(j)*(abs(u)-DATA(j,5)))/abs(u);
        else 
            Kt(j) = 0;      
        end
        Ko=Ko+Kt(j);
        end
        ff=Ko*u-F;
        %calculate du
        du = Kini^-1 *(-ff);
        residual = abs(du); % evaluate residual
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
        % Internal element force           
        for (k=1:8)
        f(k) = K(k)*u;
        end
    % Force and Dispalcement for each increment
    F1i(i) = F;
    U1(i) = u;
    DU1(i)=residual;I1(i)=i;IT1(i)=it;
    if abs(u) >= Dmax;disp('  ## Displacement reached to ultimate displacement ##');break;end
end
D1=[0;U1'];F1=[0;F1i'];
%% Linear Springs Analysis
disp('#################################################');
disp('#      Pushover Analysis of Linear Springs      #');
disp('#################################################');
% Gradually increase the applied load
for i=1:i
    % Define the applied load
        F = [P*i];
        Kini=0;
        for (j=1:8)
        K(j) = Rk1(j);      
        Kini = Kini+K(j);
        end
        
        it = 0; % initialize iteration count
        residual = 100; % initialize residual
        while (residual > tolerance)
        Ko=0;
        for (j=1:8)
        Kt(j) = Rk1(j);     
        Ko=Ko+Kt(j);
        end
        ff=Ko*u-F;
        %calculate du
        du = Kini^-1 *(-ff);
        residual = abs(du); % evaluate residual
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
        % Internal element force           
        for (k=1:8)
        f(k) = K(k)*u;
        end
% Force and Dispalcement for each increment
    F2i(i) = F;
    U2(i) = u;
    DU2(i)=residual;I2(i)=i;IT2(i)=it;
end
D2=[0;U2'];F2=[0;F2i'];
%% imaging
figure (1)
IMAGE=imread('PushoverAnalysisNonlinear8SpringsFC.jpg');
image(IMAGE);axis image;axis off;
figure(2)
p1=plot(I1,DU1,'black',I2,DU2,'--green');grid on;set(p1,'LineWidth',2);
xlabel('increment');ylabel('Residual');
legend('Nonlinear','Linear','Location','NorthEastOutside');
title('Residual-Increment diagram','color','b');
figure(3)
p1=plot(I1,IT1,'black',I2,IT2,'--green');grid on;set(p1,'LineWidth',2);
xlabel('increment');ylabel('Iteration');
legend('Nonlinear','Linear','Location','NorthEastOutside');
title('Iteration-Increment diagram','color','b');
figure(4)
p1=plot(D1,F1,'black',D2,F2,'--r');grid on;set(p1,'LineWidth',2);
legend('Nonlinear','Linear','Location','NorthEastOutside');
xlabel('Displacement');ylabel('Force');
title('Force-Displacement Diagram of Pushover Analysis Linear and Nonlinear Springs ','color','b');

