%***********************************************************%
%                 >> IN THE NAME OF ALLAH <<                %
% LINK – PLASTIC KINEMATIC LINK                             %
% Checking the analysis by Force                            %
% SAP2000:Example 6-009                                     %
%-----------------------------------------------------------%
%     This program is written by salar delavar ghashghaei   %  
%          E-mail:salar.d.ghashghaei@gmail.com              %
%***********************************************************%
clear all;close all;clc
% Define Parameters in unit: mm,kN
P1 = 0.0; % [kN]
P2 = 0.5; % [kN]
lanX = 0.0;lanY = 1.0;
Dmax = 304.8; % [mm] Max displacement
m = 800; % number of calculation
itermax = 500;% maximum number of iterations
tolerance = 1e-12; % specified tolerance for convergence
u = zeros(1,1);% initial guess value
%%% monitor cpu time
starttime = cputime;
% Element stifness first-order coefficient
%% Link Properties
d1=50.8; % [mm] Yield displacement
F1=222.4111; % [kN] Yield force
d2=152.4; % [mm]
F2=311.3755; % [kN] 
d3=406.4; % [mm] Ultimate displacement
F3=355.8578;  % [kN] Ultimate force
Rk1=(F1-0)/(d1-0);
Rk2=(F2-F1)/(d2-d1);
Rk3=(F3-F2)/(d3-d2);
%% First-order Nonlinear Analysis
disp('#################################################');
disp('#         Plastic Kinematic Link Analysis       #');
disp('#################################################');
% Gradually increase the applied load
for i=1:m
    % Define the applied load
        F = [P1;P2*i];
        if and(abs(F(2))> 0,abs(F(2))<= F1)
            K = Rk1;% link stiffness [DOF (6)]
        elseif and(abs(F(2))> F1,abs(F(2))<= F2)
            K = (F1+Rk2*(abs(u)-d1))/abs(u);% link stiffness [DOF (6)]
        elseif and(abs(F(2))> F2,abs(F(2))<= F3)
            K = (F2+Rk3*(abs(u)-d2))/abs(u);% link stiffness [DOF (6)]
        end
     Kini = [K*lanY^2];
        it = 0; % initialize iteration count
        residual = 100; % initialize residual
      while (residual > tolerance)
        % assemble global K matrix
        if and(abs(F(2))> 0,abs(F(2))<= F1)
            K = Rk1;% link stiffness [DOF (6)]
        elseif and(abs(F(2))> F1,abs(F(2))<= F2)
            K = (F1+Rk2*(abs(u)-d1))/abs(u);% link stiffness [DOF (6)]
        elseif and(abs(F(2))> F2,abs(F(2))<= F3)
            K = (F2+Rk3*(abs(u)-d2))/abs(u);% link stiffness [DOF (6)]
        end
    Ko = [K*lanY^2];
        f=Ko*u-F(2);
        %calculate du1 & du2
        du = Kini^-1 *(-f);
        %Calculate the residual (internal-external force deviation)
        residual = max(abs(du)); % evaluate residual
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
        % Stiffness Matrix for each element
        Kele1 = K*[-lanX -lanY lanX lanY]; 
        Fele1 = Kele1*[0;0;0;u];   
% Force and Dispalcement for each increment
    F1i(i) = F(2);
    U1(i) = u;
    DU(i)=residual;I(i)=i;IT(i)=it;
    INT_N_f1(i) = roundn(Fele1(1),-3);% Interal force of base Shear
    if abs(u) >= Dmax;disp('  ## Link reached to ultimate displacement ##');break;end
end
D1=[0;U1'];F1=[0;F1i'];
%% SAP2000 analysis report [mm,kN]
D2=[0
3.048
6.096
9.144
12.192
15.24
18.288
21.336
24.384
27.432
30.48
33.528
36.576
39.624
42.672
45.72
48.768
51.816
54.864
57.912
60.96
64.008
67.056
70.104
73.152
76.2
79.248
82.296
85.344
88.392
91.44
94.488
97.536
100.584
103.632
106.68
109.728
112.776
115.824
118.872
121.92
124.968
128.016
131.064
134.112
137.16
140.208
143.256
146.304
149.352
152.4
155.448
158.496
161.544
164.592
167.64
170.688
173.736
176.784
179.832
182.88
185.928
188.976
192.024
195.072
198.12
201.168
204.216
207.264
210.312
213.36
216.408
219.456
222.504
225.552
228.6
231.648
234.696
237.744
240.792
243.84
246.888
249.936
252.984
256.032
259.08
262.128
265.176
268.224
271.272
274.32
277.368
280.416
283.464
286.512
289.56
292.608
295.656
298.704
301.752
304.8];
F2=[0
13.345
26.689
40.034
53.379
66.723
80.068
93.413
106.757
120.102
133.447
146.791
160.136
173.481
186.825
200.17
213.515
223.301
225.97
228.639
231.308
233.976
236.645
239.314
241.983
244.652
247.321
249.99
252.659
255.328
257.997
260.666
263.335
266.004
268.673
271.342
274.01
276.679
279.348
282.017
284.686
287.355
290.024
292.693
295.362
298.031
300.7
303.369
306.038
308.707
311.376
311.909
312.443
312.977
313.511
314.044
314.578
315.112
315.646
316.18
316.713
317.247
317.781
318.315
318.849
319.382
319.916
320.45
320.984
321.517
322.051
322.585
323.119
323.653
324.186
324.72
325.254
325.788
326.322
326.855
327.389
327.923
328.457
328.99
329.524
330.058
330.592
331.126
331.659
332.193
332.727
333.261
333.795
334.328
334.862
335.396
335.93
336.464
336.997
337.531
338.065];
%% imaging
figure (1)
IMAGE=imread('PlasticKinematicLinkFC.jpg');
image(IMAGE);axis image;axis off;
figure(2)
p1=plot(I,DU,'black');grid on;set(p1,'LineWidth',2);
xlabel('increment');ylabel('Residual');
title('Residual-Increment diagram','color','b');
figure(3)
p1=plot(I,IT,'black');grid on;set(p1,'LineWidth',2);
xlabel('increment');ylabel('Iteration');
title('Iteration-Increment diagram','color','b');
figure(4)
p1=plot(D1,F1,D2,F2,'r--');grid on;set(p1,'LineWidth',3);
legend('MATLAB','SAP2000','Location','NorthEastOutside');
xlabel('Displacement (mm)');ylabel('Force (kN)');
title('Force-Displacement Diagram of Plastic Kinematic Link ','color','b');

