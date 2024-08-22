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
Force = [20 30 32 35];
Displacement = [5 35 50 80];
Coff(1) = .5;Coff(2) = .8;Coff(3) = 1;Coff(4)= .8;Coff(5) = .5; % all coefficients must less and equal 1
Dmax = max(Displacement); % [mm] Max displacement
Rk1=(Force(1)-0)/(Displacement(1)-0);
Rk2=(Force(2)-Force(1))/(Displacement(2)-Displacement(1));
Rk3=(Force(3)-Force(2))/(Displacement(3)-Displacement(2));
Rk4=(Force(4)-Force(3))/(Displacement(4)-Displacement(3));
f=zeros(5,1);
%% Nonlinear Springs Analysis
disp('#################################################');
disp('#    Pushover Analysis of Nonlinear Springs     #');
disp('#################################################');
% Gradually increase the applied load
for i=1:m
    % Define the applied load
        F = [P*i];
        Kini=0;
        for (j=1:5)
        if and(abs(f(j)) >= 0,abs(f(j)) <= Force(1))
            K(j) = Coff(j)*Rk1;
        elseif and(abs(f(j)) > Force(1),abs(f(j))<= Force(2))
            K(j) = Coff(j)*(Force(1)+Rk2*(abs(u)-Displacement(1)))/abs(u);
        elseif and(abs(f(j)) > Force(2),abs(f(j))<= Force(3))
            K(j) = Coff(j)*(Force(2)+Rk3*(abs(u)-Displacement(2)))/abs(u);
        elseif and(abs(f(j)) > Force(3),abs(f(j)) <= Force(4))
            K(j) = Coff(j)*(Force(3)+Rk4*(abs(u)-Displacement(3)))/abs(u);
        else 
            K(j) = 0;    
            
        end
        Kini = Kini+K(j);
        end
        
        it = 0; % initialize iteration count
        residual = 100; % initialize residual
        while (residual > tolerance)
        Ko=0;
        for (j=1:5)
        if and(abs(f(j)) >= 0,abs(f(j)) <= Force(1))
            Kt(j) = Coff(j)*Rk1;
        elseif and(abs(f(j)) > Force(1),abs(f(j))<= Force(2))
            Kt(j) = Coff(j)*(Force(1)+Rk2*(abs(u)-Displacement(1)))/abs(u);
        elseif and(abs(f(j)) > Force(2),abs(f(j))<= Force(3))
            Kt(j) = Coff(j)*(Force(2)+Rk3*(abs(u)-Displacement(2)))/abs(u);
        elseif and(abs(f(j)) > Force(3),abs(f(j)) <= Force(4))
            Kt(j) = Coff(j)*(Force(3)+Rk4*(abs(u)-Displacement(3)))/abs(u);
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
        for (k=1:5)
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
        for (j=1:5)
        K(j) = Coff(j)*Rk1;      
        Kini = Kini+K(j);
        end
        
        it = 0; % initialize iteration count
        residual = 100; % initialize residual
        while (residual > tolerance)
        Ko=0;
        for (j=1:5)
        Kt(j) = Coff(j)*Rk1;     
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
        for (k=1:5)
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
IMAGE=imread('PushoverAnalysisNonlinearSpringsFC.jpg');
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

