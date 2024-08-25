%***********************************************************%
%                  >> IN THE NAME OF GOD <<                 %
%      Moment-Curvature Analysis of 4 steel sections        %
%-----------------------------------------------------------%
%     This program is written by salar delavar ghashghaei   %  
%          E-mail:salar.d.ghashghaei@gmail.com              %
%-----------------------------------------------------------%
%Unit: Newton-Milimeter                                     %
%Given: Steel section properties ,                          %
%Calculate: Moment-Curavture                                %
% Newton-Raphson Method : Tangent procedure                 %
%***********************************************************%
%        _     __                 __                        %   
%        ^    |  |               |  |                       %
%        |    |  |               |  |                       %
%             |  |               |  |                       %
%             |  |_______________|  |  +                    %
%       bf    |   _______________   |  tw                   %
%             |  |               |  |  +                    %
%        |    |  |               |  |                       %
%        v    |  |               |  |                       %
%        -     ---               ---                        %
%             |<-   tf1+tf2+ hw   ->|                       % 
%    X                                                      %
%    ^                                                      %
%    |             (Moment - Curvature along X axis)        %
%    |                                                      %
%    +----> Y                                               %
%***********************************************************%
clear all;close all;clc;
%% Section Properties
tf1=9.2;% [mm] I section thickness on Top flange
bf1=110;% [mm] I section width on Top flange
tw=5.9;% [mm] I section thickness of Web
hw=201.6;% [mm] Height of web
tf2=9.2;% [mm] I section thickness on Bottom flange
bf2=110;% [mm] I section width on Bottom flange
ptf1=10;% [mm] Plate section thickness on Top flange [section 2]
pbf1=90;% [mm] Plate section width on Top flange [section 2]
ptf2=10;% [mm] Plate section thickness on Bottom flange [section 2]
pbf2=90;% [mm] Plate section width on Bottom flange [section 2]
ptf41=10;% [mm] Plate section thickness on Top flange [section 4]
pbf41=150;% [mm] Plate section width on Top flange [section 4]
ptf42=10;% [mm] Plate section thickness on Bottom flange [section 4]
pbf42=150;% [mm] Plate section width on Bottom flange [section 4]
%% Steel Section Properties
fy =240;% [N/mm^2] Yield strength of steel section
Es =2e5;% [N/mm^2] Modulus of elasticity of steel section
fu=1.5*fy;% Ultimate steel stress
ey=fy/Es;% Yield steel strain
esh=0.025;% Strain at steel strain-hardening
esu=0.35;% Ultimate steel strain
Esh=(fu-fy)/(esu-esh);
N=1000;% Number of steel section Fiber 
itermax = 4000;% maximum number of iterations
tolerance = 1e-12; % specified tolerance for convergence
%%% monitor cpu time
starttime = cputime;
%% ------------------ Newton Method Procedure - SECTION 1 ------------------------%
disp('#############');
disp('# SECTION 1 #');
disp('#############');
h=tf1+tf2+hw;% [mm] Height of Section
x=.25*h;% initial guess of Neuteral axis
N1=N/10;% Number of steel section top flange Fiber  
N2=8*N/10;% Number of steel section  web Fiber
N3=N/10;% Number of steel section bottom flange Fiber 
R1=(1/N1);R2=(1/N2);R3=(1/N3);
for k=1:N1;c1(k)=(.5*R1+R1*(k-1))*tf1;end % distance of each top flange Fiber 
for k=1:N2;c2(k)=tf1+(.5*R2+R2*(k-1))*hw;end % distance of each web Fiber
for k=1:N3;c3(k)=tf1+hw+(.5*R3+R3*(k-1))*tf2;end % distance of each bottom flange Fiber
c01=[c1 c2 c3];
ES01=[.2*ey .4*ey .6*ey .8*ey ey .2*esh .4*esh .6*esh .8*esh esh .2*esu .4*esu .6*esu .8*esu esu];q01=size(ES01,2);
% ES01=ey:.001:esu;q01=size(ES01,2);
for j=1:q01;
    eS=ES01(j);
    it = 0; % initialize iteration count
    residual = 100; % initialize residual
    while (residual > tolerance)
        for z=1:N;% in this step: steel section force for each fiber is calculated
        es=eS*(x-c01(z))/x;            
        if and(z>=1,z<=roundn(N/10,0)) % TOP SECTION
    %---------------- Fs -------------%
       if and(es>0,es<=ey)
        fs=Es*es;
        fstan=(Es*eS*c01(z))/(x)^2;
       elseif and(es<0,es>=(-ey))
        fs=Es*es;
        fstan=(Es*eS*c01(z))/(x)^2;
       elseif  and(es>ey,es<=esh)
        fs=fy;
        fstan=0;
       elseif and(es<(-ey),es>=(-esh))
        fs=-fy;
        fstan=0;
        elseif  and(es>esh,es<=esu)
        fs=fu-(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c01(z)*(fu-fy)*(((eS*c01(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif and(es<(-esh),es>=(-esu))
        fs=-fu+(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c01(z)*(fu-fy)*(((eS*c01(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif or(es>esu,es<(-esu))
        fs=0;
        fstan=0;
       end
     Fs(z)=bf1*R1*tf1*fs;% force on top flange
     Fstan(z)=bf1*R1*tf1*fstan;% tangent steel force
        elseif and(z>roundn(N/10,0),z<=roundn(9*N/10,0)) % MIDDLE SECTION
       if and(es>0,es<=ey)
        fs=Es*es;
        fstan=(Es*eS*c01(z))/(x)^2;
       elseif and(es<0,es>=(-ey))
        fs=Es*es;
        fstan=(Es*eS*c01(z))/(x)^2;
       elseif  and(es>ey,es<=esh)
        fs=fy;
        fstan=0;
       elseif and(es<(-ey),es>=(-esh))
        fs=-fy;
        fstan=0;
        elseif  and(es>esh,es<=esu)
        fs=fu-(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c01(z)*(fu-fy)*(((eS*c01(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif and(es<(-esh),es>=(-esu))
        fs=-fu+(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c01(z)*(fu-fy)*(((eS*c01(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif or(es>esu,es<(-esu))
        fs=0;
        fstan=0;
       end
     Fs(z)=tw*R2*hw*fs;% force on web
     Fstan(z)=tw*R2*hw*fstan;% tangent steel force
     elseif and(z>roundn(9*N/10,0),z<=N) % BELOW SECTION
       if and(es>0,es<=ey)
        fs=Es*es;
        fstan=(Es*eS*c01(z))/(x)^2;
       elseif and(es<0,es>=(-ey))
        fs=Es*es;
        fstan=(Es*eS*c01(z))/(x)^2;
       elseif  and(es>ey,es<=esh)
        fs=fy;
        fstan=0;
       elseif and(es<(-ey),es>=(-esh))
        fs=-fy;
        fstan=0;
        elseif  and(es>esh,es<=esu)
        fs=fu-(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c01(z)*(fu-fy)*(((eS*c01(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif and(es<(-esh),es>=(-esu))
        fs=-fu+(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c01(z)*(fu-fy)*(((eS*c01(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif or(es>esu,es<(-esu))
        fs=0;
        fstan=0;
       end
      Fs(z)=bf2*R3*tf2*fs;% force on bottom flange
      Fstan(z)=bf2*R3*tf2*fstan;% tangent steel force
        end
     Ss01(z)=es;SS01(j,z)=Ss01(z);%steel sction Fiber Strain
     CFS01(j,z)=fs;%steel sction Fiber Stress
        end
    %----------------------------------%
    FsTOTAL=sum(Fs);A=FsTOTAL;
    FsTOTAL_tan=sum(Fstan);A_tan=FsTOTAL_tan;
    dx = A_tan^-1 *(-A);
    residual = abs(dx); % evaluate residual
        it = it + 1; % increment iteration count
        x = x+dx; % update x
        if it == itermax % stop the the analysis of this step please of Convergence
          fprintf('(-)Increment %1.0f : trail iteration reached to Ultimate %1.0f - strain: %1.6f - error: [%1.2f]\n',j,it,eS,A)
             disp('    ## The solution for this step is not converged. Please check your model ##') 
            break
        end
    end
    if it == itermax;break;end % stop the analysis at all because last Convergence     
     for k=1:N;cc01(k)=x-c01(k);end %distance of each steel fiber from Neuteral axis
        if it < itermax% iteration control
        fprintf('(+)Increment %1.0f : It is converged in %1.0f iterations - strain: %1.5f - x: %1.2f - Phi: %1.5f - Moment: %1.2f\n',j,it,eS,x,(eS/x)*1000,(Fs*cc01')*10^-6)
        end
 if and(eS>0,eS<ey);ESS01(j)=eS;TUSS01(j)=Es*eS;elseif and(eS>=ey,eS<esh);ESS01(j)=eS;TUSS01(j)=fy;elseif and(eS>=esh,eS<esu);ESS01(j)=eS;TUSS01(j)=fu-(fu-fy)*((esu-abs(eS))/(esu-esh))^2;end% Top Steel section compression stress        
% Calculate Moment and Curavture
Cur01(j)=(eS/x)*1000;XX01(j)=x;AA01(j)=A;
Mom01(j)=(Fs*cc01')*10^-6;
CUR01(j)=Cur01(j);I01(j)=j;IT01(j)=it;DX01(j)=dx;
end
Cur01=[0 Cur01];Mom01=[0 Mom01];
s01=size(Cur01,2);for i=1:s01-1;EI01(i)=(Mom01(i+1)-Mom01(i))/(Cur01(i+1)-Cur01(i));end % Flextural Rigidity
if roundn(eS,-5) == esu;fprintf('\n      ## Strain Reached to Ultimate Strain: %1.4f ## \n\n',eS);end
%% Steel Section bilinear fitting - SECTION 1
SIZE01=size(Mom01,2);
for i=1:SIZE01-1;
    hh01(i) = Cur01(i+1)-Cur01(i);
    Aa01(i)=(Mom01(i)+Mom01(i+1))*0.5*hh01(i);
end
Area01=sum(Aa01);k001 =Mom01(2)/Cur01(2);
fiy01 = (Mom01(i+1)*max(Cur01)*0.5-Area01)/(Mom01(i+1)*0.5 - k001*max(Cur01)*0.5);
My01 = k001*fiy01;
X01 = [0 fiy01 max(Cur01)];Y01 = [0 My01 Mom01(i+1)];
disp('+=============================+');
disp('=   Section 1 curve fitted    =');
disp('     Curvature    Moment       ');
disp('       (1/m)      (kN.m)       ');
disp('-------------------------------');
disp([X01' Y01']);
disp('+=============================+');
%% ------------------ Newton Method Procedure - SECTION 2 ------------------------%
disp('#############');
disp('# SECTION 2 #');
disp('#############');
clear c1 c2 c3
h=ptf1+ptf2+tf1+tf2+hw;% [mm] Height of Section
x=.25*h;% initial guess of Neuteral axis
N1=N/10;% Number of I steel section top plate flange Fiber
N2=N/10;% Number of steel section top flange Fiber 
N3=6*N/10;% Number of I steel section web Fiber
N4=N/10;% Number of I steel section bottom flange Fiber
N5=N/10;% Number of steel section bottom Plate flange Fiber 
R1=(1/N1);R2=(1/N2);R3=(1/N3);R4=(1/N4);R5=(1/N5);
for k=1:N1;c1(k)=(.5*R1+R1*(k-1))*ptf1;end % distance of each top Plate flange Fiber 
for k=1:N2;c2(k)=ptf1+(.5*R2+R2*(k-1))*tf1;end % distance of each top I flange Fiber
for k=1:N3;c3(k)=ptf1+tf1+(.5*R3+R3*(k-1))*hw;end % distance of each web Fiber
for k=1:N4;c4(k)=ptf1+tf1+hw+(.5*R4+R4*(k-1))*tf2;end % distance of each bottom I flange Fiber
for k=1:N5;c5(k)=ptf1+tf1+hw+tf2+(.5*R5+R5*(k-1))*ptf2;end % distance of each bottom Plate flange Fiber
c02=[c1 c2 c3 c4 c5];
ES02=[.2*ey .4*ey .6*ey .8*ey ey .2*esh .4*esh .6*esh .8*esh esh .2*esu .4*esu .6*esu .8*esu esu];q02=size(ES02,2);
% ES02=ey:.001:esu;q02=size(ES02,2);
for j=1:q02;
    eS=ES02(j);
    it = 0; % initialize iteration count
    residual = 100; % initialize residual
    while (residual > tolerance)
        for z=1:N;% in this step: steel section force for each fiber is calculated
        es=eS*(x-c02(z))/x;            
        if and(z>=1,z<=roundn(N/10,0)) % Top Plate of flange
    %---------------- Fs -------------%
       if and(es>0,es<=ey)
        fs=Es*es;
        fstan=(Es*eS*c02(z))/(x)^2;
       elseif and(es<0,es>=(-ey))
        fs=Es*es;
        fstan=(Es*eS*c02(z))/(x)^2;
       elseif  and(es>ey,es<=esh)
        fs=fy;
        fstan=0;
       elseif and(es<(-ey),es>=(-esh))
        fs=-fy;
        fstan=0;
        elseif  and(es>esh,es<=esu)
        fs=fu-(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c02(z)*(fu-fy)*(((eS*c02(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif and(es<(-esh),es>=(-esu))
        fs=-fu+(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c02(z)*(fu-fy)*(((eS*c02(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif or(es>esu,es<(-esu))
        fs=0;
        fstan=0;
       end
     Fs(z)=pbf1*R1*ptf1*fs;% force on top flange
     Fstan(z)=pbf1*R1*ptf1*fstan;% tangent steel force
        elseif and(z>roundn(N/10,0),z<=roundn(2*N/10,0)) % Top flange of I
       if and(es>0,es<=ey)
        fs=Es*es;
        fstan=(Es*eS*c02(z))/(x)^2;
       elseif and(es<0,es>=(-ey))
        fs=Es*es;
        fstan=(Es*eS*c02(z))/(x)^2;
       elseif  and(es>ey,es<=esh)
        fs=fy;
        fstan=0;
       elseif and(es<(-ey),es>=(-esh))
        fs=-fy;
        fstan=0;
        elseif  and(es>esh,es<=esu)
        fs=fu-(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c02(z)*(fu-fy)*(((eS*c02(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif and(es<(-esh),es>=(-esu))
        fs=-fu+(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c02(z)*(fu-fy)*(((eS*c02(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif or(es>esu,es<(-esu))
        fs=0;
        fstan=0;
       end
     Fs(z)=bf1*R2*tf1*fs;% force on top flange
     Fstan(z)=bf1*R2*tf1*fstan;% tangent steel force
        elseif and(z>roundn(2*N/10,0),z<=roundn(8*N/10,0)) % Middle web
       if and(es>0,es<=ey)
        fs=Es*es;
        fstan=(Es*eS*c02(z))/(x)^2;
       elseif and(es<0,es>=(-ey))
        fs=Es*es;
        fstan=(Es*eS*c02(z))/(x)^2;
       elseif  and(es>ey,es<=esh)
        fs=fy;
        fstan=0;
       elseif and(es<(-ey),es>=(-esh))
        fs=-fy;
        fstan=0;
        elseif  and(es>esh,es<=esu)
        fs=fu-(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c02(z)*(fu-fy)*(((eS*c02(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif and(es<(-esh),es>=(-esu))
        fs=-fu+(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c02(z)*(fu-fy)*(((eS*c02(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif or(es>esu,es<(-esu))
        fs=0;
        fstan=0;
       end
     Fs(z)=tw*R3*hw*fs;% force on web
     Fstan(z)=tw*R3*hw*fstan;% tangent steel force
     elseif and(z>roundn(8*N/10,0),z<=roundn(9*N/10,0)) % Below flange of I
       if and(es>0,es<=ey)
        fs=Es*es;
        fstan=(Es*eS*c02(z))/(x)^2;
       elseif and(es<0,es>=(-ey))
        fs=Es*es;
        fstan=(Es*eS*c02(z))/(x)^2;
       elseif  and(es>ey,es<=esh)
        fs=fy;
        fstan=0;
       elseif and(es<(-ey),es>=(-esh))
        fs=-fy;
        fstan=0;
        elseif  and(es>esh,es<=esu)
        fs=fu-(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c02(z)*(fu-fy)*(((eS*c02(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif and(es<(-esh),es>=(-esu))
        fs=-fu+(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c02(z)*(fu-fy)*(((eS*c02(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif or(es>esu,es<(-esu))
        fs=0;
        fstan=0;
       end
      Fs(z)=bf2*R4*tf2*fs;% force on bottom flange
      Fstan(z)=bf2*R4*tf2*fstan;% tangent steel force
      elseif and(z>roundn(9*N/10,0),z<=N) % below Plate of flange
       if and(es>0,es<=ey)
        fs=Es*es;
        fstan=(Es*eS*c02(z))/(x)^2;
       elseif and(es<0,es>=(-ey))
        fs=Es*es;
        fstan=(Es*eS*c02(z))/(x)^2;
       elseif  and(es>ey,es<=esh)
        fs=fy;
        fstan=0;
       elseif and(es<(-ey),es>=(-esh))
        fs=-fy;
        fstan=0;
        elseif  and(es>esh,es<=esu)
        fs=fu-(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c02(z)*(fu-fy)*(((eS*c02(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif and(es<(-esh),es>=(-esu))
        fs=-fu+(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c02(z)*(fu-fy)*(((eS*c02(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif or(es>esu,es<(-esu))
        fs=0;
        fstan=0;
       end
      Fs(z)=pbf2*R5*ptf2*fs;% force on bottom flange
      Fstan(z)=pbf2*R5*ptf2*fstan;% tangent steel force
        end
     Ss02(z)=es;SS02(j,z)=Ss02(z);%steel sction Fiber Strain
     CFS02(j,z)=fs;%steel sction Fiber Stress
        end
    %----------------------------------%
    FsTOTAL=sum(Fs);A=FsTOTAL;
    FsTOTAL_tan=sum(Fstan);A_tan=FsTOTAL_tan;
    dx = A_tan^-1 *(-A);
    residual = abs(dx); % evaluate residual
        it = it + 1; % increment iteration count
        x = x+dx; % update x
        if it == itermax % stop the the analysis of this step please of Convergence
          fprintf('(-)Increment %1.0f : trail iteration reached to Ultimate %1.0f - strain: %1.6f - error: [%1.2f]\n',j,it,eS,A)
             disp('    ## The solution for this step is not converged. Please check your model ##') 
            break
        end
    end
    if it == itermax;break;end % stop the analysis at all because last Convergence     
     for k=1:N;cc02(k)=x-c02(k);end %distance of each steel fiber from Neuteral axis
        if it < itermax% iteration control
        fprintf('(+)Increment %1.0f : It is converged in %1.0f iterations - strain: %1.5f - x: %1.2f - Phi: %1.5f - Moment: %1.2f\n',j,it,eS,x,(eS/x)*1000,(Fs*cc02')*10^-6)
        end
 if and(eS>0,eS<ey);ESS02(j)=eS;TUSS02(j)=Es*eS;elseif and(eS>=ey,eS<esh);ESS02(j)=eS;TUSS02(j)=fy;elseif and(eS>=esh,eS<esu);ESS02(j)=eS;TUSS02(j)=fu-(fu-fy)*((esu-abs(eS))/(esu-esh))^2;end% Top Steel section compression stress        
% Calculate Moment and Curavture
Cur02(j)=(eS/x)*1000;XX02(j)=x;AA02(j)=A;
Mom02(j)=(Fs*cc02')*10^-6;
CUR02(j)=Cur02(j);I02(j)=j;IT02(j)=it;DX02(j)=dx;
end
if roundn(eS,-5) == esu;fprintf('\n      ## Strain Reached to Ultimate Strain: %1.4f ## \n\n',eS);end
Cur02=[0 Cur02];Mom02=[0 Mom02];
s02=size(Cur02,2);for i=1:s02-1;EI02(i)=(Mom02(i+1)-Mom02(i))/(Cur02(i+1)-Cur02(i));end % Flextural Rigidity
%% Steel Section bilinear fitting - SECTION 2
SIZE02=size(Mom02,2);
for i=1:SIZE02-1;
    hh02(i) = Cur02(i+1)-Cur02(i);
    Aa02(i)=(Mom02(i)+Mom02(i+1))*0.5*hh02(i);
end
Area02=sum(Aa02);k002 =Mom02(2)/Cur02(2);
fiy02 = (Mom02(i+1)*max(Cur02)*0.5-Area02)/(Mom02(i+1)*0.5 - k002*max(Cur02)*0.5);
My02 = k002*fiy02;
X02 = [0 fiy02 max(Cur02)];Y02 = [0 My02 Mom02(i+1)];
disp('+=============================+');
disp('=   Section 2 curve fitted    =');
disp('     Curvature    Moment       ');
disp('       (1/m)      (kN.m)       ');
disp('-------------------------------');
disp([X02' Y02']);
disp('+=============================+');
%% ------------------ Newton Method Procedure - SECTION 3 ------------------------%
disp('#############');
disp('# SECTION 3 #');
disp('#############');
clear c1 c2 c3 c4 c5
h=tf1+tf2+hw;% [mm] Height of Section
x=.25*h;% initial guess of Neuteral axis
N1=N/10;% Number of steel section top flange Fiber  
N2=8*N/10;% Number of steel section  web Fiber
N3=N/10;% Number of steel section bottom flange Fiber 
R1=(1/N1);R2=(1/N2);R3=(1/N3);
for k=1:N1;c1(k)=(.5*R1+R1*(k-1))*tf1;end % distance of each top flange Fiber 
for k=1:N2;c2(k)=tf1+(.5*R2+R2*(k-1))*hw;end % distance of each web Fiber
for k=1:N3;c3(k)=tf1+hw+(.5*R3+R3*(k-1))*tf2;end % distance of each bottom flange Fiber
c03=[c1 c2 c3];
ES03=[.2*ey .4*ey .6*ey .8*ey ey .2*esh .4*esh .6*esh .8*esh esh .2*esu .4*esu .6*esu .8*esu esu];q03=size(ES03,2);
% ES03=ey:.001:esu;q03=size(ES03,2);
for j=1:q03;
    eS=ES03(j);
    it = 0; % initialize iteration count
    residual = 100; % initialize residual
    while (residual > tolerance)
        for z=1:N;% in this step: steel section force for each fiber is calculated
        es=eS*(x-c03(z))/x;            
        if and(z>=1,z<=roundn(N/10,0)) % TOP SECTION
    %---------------- Fs -------------%
       if and(es>0,es<=ey)
        fs=Es*es;
        fstan=(Es*eS*c03(z))/(x)^2;
       elseif and(es<0,es>=(-ey))
        fs=Es*es;
        fstan=(Es*eS*c03(z))/(x)^2;
       elseif  and(es>ey,es<=esh)
        fs=fy;
        fstan=0;
       elseif and(es<(-ey),es>=(-esh))
        fs=-fy;
        fstan=0;
        elseif  and(es>esh,es<=esu)
        fs=fu-(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c03(z)*(fu-fy)*(((eS*c03(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif and(es<(-esh),es>=(-esu))
        fs=-fu+(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c03(z)*(fu-fy)*(((eS*c03(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif or(es>esu,es<(-esu))
        fs=0;
        fstan=0;
       end
     Fs(z)=2*bf1*R1*tf1*fs;% force on top flange
     Fstan(z)=2*bf1*R1*tf1*fstan;% tangent steel force
        elseif and(z>roundn(N/10,0),z<=roundn(9*N/10,0)) % MIDDLE SECTION
       if and(es>0,es<=ey)
        fs=Es*es;
        fstan=(Es*eS*c03(z))/(x)^2;
       elseif and(es<0,es>=(-ey))
        fs=Es*es;
        fstan=(Es*eS*c03(z))/(x)^2;
       elseif  and(es>ey,es<=esh)
        fs=fy;
        fstan=0;
       elseif and(es<(-ey),es>=(-esh))
        fs=-fy;
        fstan=0;
        elseif  and(es>esh,es<=esu)
        fs=fu-(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c03(z)*(fu-fy)*(((eS*c03(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif and(es<(-esh),es>=(-esu))
        fs=-fu+(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c03(z)*(fu-fy)*(((eS*c03(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif or(es>esu,es<(-esu))
        fs=0;
        fstan=0;
       end
     Fs(z)=2*tw*R2*hw*fs;% force on web
     Fstan(z)=2*tw*R2*hw*fstan;% tangent steel force
     elseif and(z>roundn(9*N/10,0),z<=N) % BELOW SECTION
       if and(es>0,es<=ey)
        fs=Es*es;
        fstan=(Es*eS*c03(z))/(x)^2;
       elseif and(es<0,es>=(-ey))
        fs=Es*es;
        fstan=(Es*eS*c03(z))/(x)^2;
       elseif  and(es>ey,es<=esh)
        fs=fy;
        fstan=0;
       elseif and(es<(-ey),es>=(-esh))
        fs=-fy;
        fstan=0;
        elseif  and(es>esh,es<=esu)
        fs=fu-(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c03(z)*(fu-fy)*(((eS*c03(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif and(es<(-esh),es>=(-esu))
        fs=-fu+(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c03(z)*(fu-fy)*(((eS*c03(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif or(es>esu,es<(-esu))
        fs=0;
        fstan=0;
       end
      Fs(z)=2*bf2*R3*tf2*fs;% force on bottom flange
      Fstan(z)=2*bf2*R3*tf2*fstan;% tangent steel force
        end
     Ss03(z)=es;SS03(j,z)=Ss03(z);%steel sction Fiber Strain
     CFS03(j,z)=fs;%steel sction Fiber Stress
        end
    %----------------------------------%
    FsTOTAL=sum(Fs);A=FsTOTAL;
    FsTOTAL_tan=sum(Fstan);A_tan=FsTOTAL_tan;
    dx = A_tan^-1 *(-A);
    residual = abs(dx); % evaluate residual
        it = it + 1; % increment iteration count
        x = x+dx; % update x
        if it == itermax % stop the the analysis of this step please of Convergence
          fprintf('(-)Increment %1.0f : trail iteration reached to Ultimate %1.0f - strain: %1.6f - error: [%1.2f]\n',j,it,eS,A)
             disp('    ## The solution for this step is not converged. Please check your model ##') 
            break
        end
    end
    if it == itermax;break;end % stop the analysis at all because last Convergence     
     for k=1:N;cc03(k)=x-c03(k);end %distance of each steel fiber from Neuteral axis
        if it < itermax% iteration control
        fprintf('(+)Increment %1.0f : It is converged in %1.0f iterations - strain: %1.5f - x: %1.2f - Phi: %1.5f - Moment: %1.2f\n',j,it,eS,x,(eS/x)*1000,(Fs*cc03')*10^-6)
        end
 if and(eS>0,eS<ey);ESS03(j)=eS;TUSS03(j)=Es*eS;elseif and(eS>=ey,eS<esh);ESS03(j)=eS;TUSS03(j)=fy;elseif and(eS>=esh,eS<esu);ESS03(j)=eS;TUSS03(j)=fu-(fu-fy)*((esu-abs(eS))/(esu-esh))^2;end% Top Steel section compression stress        
% Calculate Moment and Curavture
Cur03(j)=(eS/x)*1000;XX03(j)=x;AA03(j)=A;
Mom03(j)=(Fs*cc03')*10^-6;
CUR03(j)=Cur03(j);I03(j)=j;IT03(j)=it;DX03(j)=dx;
end
Cur03=[0 Cur03];Mom03=[0 Mom03];
s03=size(Cur03,2);for i=1:s03-1;EI03(i)=(Mom03(i+1)-Mom03(i))/(Cur03(i+1)-Cur03(i));end % Flextural Rigidity
if roundn(eS,-5) == esu;fprintf('\n      ## Strain Reached to Ultimate Strain: %1.4f ## \n\n',eS);end
%% Steel Section bilinear fitting - SECTION 3
SIZE03=size(Mom03,2);
for i=1:SIZE03-1;
    hh03(i) = Cur03(i+1)-Cur03(i);
    Aa03(i)=(Mom03(i)+Mom03(i+1))*0.5*hh03(i);
end
Area03=sum(Aa03);k003 =Mom03(2)/Cur03(2);
fiy03 = (Mom03(i+1)*max(Cur03)*0.5-Area03)/(Mom03(i+1)*0.5 - k003*max(Cur03)*0.5);
My03 = k003*fiy03;
X03 = [0 fiy03 max(Cur03)];Y03 = [0 My03 Mom03(i+1)];
disp('+=============================+');
disp('=   Section 3 curve fitted    =');
disp('     Curvature    Moment       ');
disp('       (1/m)      (kN.m)       ');
disp('-------------------------------');
disp([X03' Y03']);
disp('+=============================+');
%% ------------------ Newton Method Procedure - SECTION 4 ------------------------%
disp('#############');
disp('# SECTION 4 #');
disp('#############');
clear c1 c2 c3
h=ptf1+ptf2+tf1+tf2+hw;% [mm] Height of Section
x=.25*h;% initial guess of Neuteral axis
N1=N/10;% Number of I steel section top plate flange Fiber
N2=N/10;% Number of steel section top flange Fiber 
N3=6*N/10;% Number of I steel section web Fiber
N4=N/10;% Number of I steel section bottom flange Fiber
N5=N/10;% Number of steel section bottom Plate flange Fiber 
R1=(1/N1);R2=(1/N2);R3=(1/N3);R4=(1/N4);R5=(1/N5);
for k=1:N1;c1(k)=(.5*R1+R1*(k-1))*ptf1;end % distance of each top Plate flange Fiber 
for k=1:N2;c2(k)=ptf1+(.5*R2+R2*(k-1))*tf1;end % distance of each top I flange Fiber
for k=1:N3;c3(k)=ptf1+tf1+(.5*R3+R3*(k-1))*hw;end % distance of each web Fiber
for k=1:N4;c4(k)=ptf1+tf1+hw+(.5*R4+R4*(k-1))*tf2;end % distance of each bottom I flange Fiber
for k=1:N5;c5(k)=ptf1+tf1+hw+tf2+(.5*R5+R5*(k-1))*ptf2;end % distance of each bottom Plate flange Fiber
c04=[c1 c2 c3 c4 c5];
ES04=[.2*ey .4*ey .6*ey .8*ey ey .2*esh .4*esh .6*esh .8*esh esh .2*esu .4*esu .6*esu .8*esu esu];q04=size(ES04,2);
% ES04=ey:.001:esu;q04=size(ES04,2);
for j=1:q04;
    eS=ES04(j);
    it = 0; % initialize iteration count
    residual = 100; % initialize residual
    while (residual > tolerance)
        for z=1:N;% in this step: steel section force for each fiber is calculated
        es=eS*(x-c04(z))/x;            
        if and(z>=1,z<=roundn(N/10,0)) % Top Plate of flange
    %---------------- Fs -------------%
       if and(es>0,es<=ey)
        fs=Es*es;
        fstan=(Es*eS*c04(z))/(x)^2;
       elseif and(es<0,es>=(-ey))
        fs=Es*es;
        fstan=(Es*eS*c04(z))/(x)^2;
       elseif  and(es>ey,es<=esh)
        fs=fy;
        fstan=0;
       elseif and(es<(-ey),es>=(-esh))
        fs=-fy;
        fstan=0;
        elseif  and(es>esh,es<=esu)
        fs=fu-(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c04(z)*(fu-fy)*(((eS*c04(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif and(es<(-esh),es>=(-esu))
        fs=-fu+(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c04(z)*(fu-fy)*(((eS*c04(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif or(es>esu,es<(-esu))
        fs=0;
        fstan=0;
       end
     Fs(z)=pbf41*R1*ptf41*fs;% force on top flange
     Fstan(z)=pbf41*R1*ptf41*fstan;% tangent steel force
        elseif and(z>roundn(N/10,0),z<=roundn(2*N/10,0)) % Top flange of I
       if and(es>0,es<=ey)
        fs=Es*es;
        fstan=(Es*eS*c04(z))/(x)^2;
       elseif and(es<0,es>=(-ey))
        fs=Es*es;
        fstan=(Es*eS*c04(z))/(x)^2;
       elseif  and(es>ey,es<=esh)
        fs=fy;
        fstan=0;
       elseif and(es<(-ey),es>=(-esh))
        fs=-fy;
        fstan=0;
        elseif  and(es>esh,es<=esu)
        fs=fu-(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c04(z)*(fu-fy)*(((eS*c04(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif and(es<(-esh),es>=(-esu))
        fs=-fu+(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c04(z)*(fu-fy)*(((eS*c04(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif or(es>esu,es<(-esu))
        fs=0;
        fstan=0;
       end
     Fs(z)=2*bf1*R2*tf1*fs;% force on top flange
     Fstan(z)=2*bf1*R2*tf1*fstan;% tangent steel force
        elseif and(z>roundn(2*N/10,0),z<=roundn(8*N/10,0)) % Middle web
       if and(es>0,es<=ey)
        fs=Es*es;
        fstan=(Es*eS*c04(z))/(x)^2;
       elseif and(es<0,es>=(-ey))
        fs=Es*es;
        fstan=(Es*eS*c04(z))/(x)^2;
       elseif  and(es>ey,es<=esh)
        fs=fy;
        fstan=0;
       elseif and(es<(-ey),es>=(-esh))
        fs=-fy;
        fstan=0;
        elseif  and(es>esh,es<=esu)
        fs=fu-(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c04(z)*(fu-fy)*(((eS*c04(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif and(es<(-esh),es>=(-esu))
        fs=-fu+(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c04(z)*(fu-fy)*(((eS*c04(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif or(es>esu,es<(-esu))
        fs=0;
        fstan=0;
       end
     Fs(z)=2*tw*R3*hw*fs;% force on web
     Fstan(z)=2*tw*R3*hw*fstan;% tangent steel force
     elseif and(z>roundn(8*N/10,0),z<=roundn(9*N/10,0)) % Below flange of I
       if and(es>0,es<=ey)
        fs=Es*es;
        fstan=(Es*eS*c04(z))/(x)^2;
       elseif and(es<0,es>=(-ey))
        fs=Es*es;
        fstan=(Es*eS*c04(z))/(x)^2;
       elseif  and(es>ey,es<=esh)
        fs=fy;
        fstan=0;
       elseif and(es<(-ey),es>=(-esh))
        fs=-fy;
        fstan=0;
        elseif  and(es>esh,es<=esu)
        fs=fu-(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c04(z)*(fu-fy)*(((eS*c04(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif and(es<(-esh),es>=(-esu))
        fs=-fu+(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c04(z)*(fu-fy)*(((eS*c04(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif or(es>esu,es<(-esu))
        fs=0;
        fstan=0;
       end
      Fs(z)=2*bf2*R4*tf2*fs;% force on bottom flange
      Fstan(z)=2*bf2*R4*tf2*fstan;% tangent steel force
      elseif and(z>roundn(9*N/10,0),z<=N) % below Plate of flange
       if and(es>0,es<=ey)
        fs=Es*es;
        fstan=(Es*eS*c04(z))/(x)^2;
       elseif and(es<0,es>=(-ey))
        fs=Es*es;
        fstan=(Es*eS*c04(z))/(x)^2;
       elseif  and(es>ey,es<=esh)
        fs=fy;
        fstan=0;
       elseif and(es<(-ey),es>=(-esh))
        fs=-fy;
        fstan=0;
        elseif  and(es>esh,es<=esu)
        fs=fu-(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c04(z)*(fu-fy)*(((eS*c04(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif and(es<(-esh),es>=(-esu))
        fs=-fu+(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c04(z)*(fu-fy)*(((eS*c04(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif or(es>esu,es<(-esu))
        fs=0;
        fstan=0;
       end
      Fs(z)=pbf42*R5*ptf42*fs;% force on bottom flange
      Fstan(z)=pbf42*R5*ptf42*fstan;% tangent steel force
        end
     Ss04(z)=es;SS04(j,z)=Ss04(z);%steel sction Fiber Strain
     CFS04(j,z)=fs;%steel sction Fiber Stress
        end
    %----------------------------------%
    FsTOTAL=sum(Fs);A=FsTOTAL;
    FsTOTAL_tan=sum(Fstan);A_tan=FsTOTAL_tan;
    dx = A_tan^-1 *(-A);
    residual = abs(dx); % evaluate residual
        it = it + 1; % increment iteration count
        x = x+dx; % update x
        if it == itermax % stop the the analysis of this step please of Convergence
          fprintf('(-)Increment %1.0f : trail iteration reached to Ultimate %1.0f - strain: %1.6f - error: [%1.2f]\n',j,it,eS,A)
             disp('    ## The solution for this step is not converged. Please check your model ##') 
            break
        end
    end
    if it == itermax;break;end % stop the analysis at all because last Convergence     
     for k=1:N;cc04(k)=x-c04(k);end %distance of each steel fiber from Neuteral axis
        if it < itermax% iteration control
        fprintf('(+)Increment %1.0f : It is converged in %1.0f iterations - strain: %1.5f - x: %1.2f - Phi: %1.5f - Moment: %1.2f\n',j,it,eS,x,(eS/x)*1000,(Fs*cc04')*10^-6)
        end
 if and(eS>0,eS<ey);ESS04(j)=eS;TUSS04(j)=Es*eS;elseif and(eS>=ey,eS<esh);ESS04(j)=eS;TUSS04(j)=fy;elseif and(eS>=esh,eS<esu);ESS04(j)=eS;TUSS04(j)=fu-(fu-fy)*((esu-abs(eS))/(esu-esh))^2;end% Top Steel section compression stress        
% Calculate Moment and Curavture
Cur04(j)=(eS/x)*1000;XX04(j)=x;AA04(j)=A;
Mom04(j)=(Fs*cc04')*10^-6;
CUR04(j)=Cur04(j);I04(j)=j;IT04(j)=it;DX04(j)=dx;
end
if roundn(eS,-5) == esu;fprintf('\n      ## Strain Reached to Ultimate Strain: %1.4f ## \n\n',eS);end
Cur04=[0 Cur04];Mom04=[0 Mom04];
s04=size(Cur04,2);for i=1:s04-1;EI04(i)=(Mom04(i+1)-Mom04(i))/(Cur04(i+1)-Cur04(i));end % Flextural Rigidity
%% Steel Section bilinear fitting - SECTION 4
SIZE04=size(Mom04,2);
for i=1:SIZE04-1;
    hh04(i) = Cur04(i+1)-Cur04(i);
    Aa04(i)=(Mom04(i)+Mom04(i+1))*0.5*hh04(i);
end
Area04=sum(Aa04);k004 =Mom04(2)/Cur04(2);
fiy04 = (Mom04(i+1)*max(Cur04)*0.5-Area04)/(Mom04(i+1)*0.5 - k004*max(Cur04)*0.5);
My04 = k004*fiy04;
X04 = [0 fiy04 max(Cur04)];Y04 = [0 My04 Mom04(i+1)];
disp('+=============================+');
disp('=   Section 2 curve fitted    =');
disp('     Curvature    Moment       ');
disp('       (1/m)      (kN.m)       ');
disp('-------------------------------');
disp([X04' Y04']);
disp('+=============================+');
%% EI and Ductility_Rito of Steel Sections
Elastic_EI01=Y01(2)/X01(2);
Plastic_EI01=(Y01(3)-Y01(2))/(X01(3)-X01(2));
Ductility_Rito01=X01(3)/X01(2);
Over_Strength_Factor01=Y01(3)/Y01(2);
Elastic_EI02=Y02(2)/X02(2);
Plastic_EI02=(Y02(3)-Y02(2))/(X02(3)-X02(2));
Ductility_Rito02=X02(3)/X02(2);
Over_Strength_Factor02=Y02(3)/Y02(2);
Elastic_EI03=Y03(2)/X03(2);
Plastic_EI03=(Y03(3)-Y03(2))/(X03(3)-X03(2));
Ductility_Rito03=X03(3)/X03(2);
Over_Strength_Factor03=Y03(3)/Y03(2);
Elastic_EI04=Y04(2)/X04(2);
Plastic_EI04=(Y04(3)-Y04(2))/(X04(3)-X04(2));
Ductility_Rito04=X04(3)/X04(2);
Over_Strength_Factor04=Y04(3)/Y04(2);
Material_Ductility_Rito=esu/esh;
fprintf('+-------------------------------------------------------+\n')
fprintf(' SECTION 1 - Elastic EI : %5.2f (kN.m^2)\n',Elastic_EI01)
fprintf(' SECTION 1 - Plastic EI : %5.2f (kN.m^2)\n',Plastic_EI01)
fprintf(' SECTION 2 - Elastic EI : %5.2f (kN.m^2)\n',Elastic_EI02)
fprintf(' SECTION 2 - Plastic EI : %5.2f (kN.m^2)\n',Plastic_EI02)
fprintf(' SECTION 3 - Elastic EI : %5.2f (kN.m^2)\n',Elastic_EI03)
fprintf(' SECTION 3 - Plastic EI : %5.2f (kN.m^2)\n',Plastic_EI03)
fprintf(' SECTION 4 - Elastic EI : %5.2f (kN.m^2)\n',Elastic_EI04)
fprintf(' SECTION 4 - Plastic EI : %5.2f (kN.m^2)\n',Plastic_EI04)
fprintf(' Steel Material Ductility Rito : %5.2f\n',Material_Ductility_Rito)
fprintf(' SECTION 1 - Steel Section Ductility Rito : %5.2f\n',Ductility_Rito01)
fprintf(' SECTION 1 - Steel Section Over Strength Factor : %5.2f\n',Over_Strength_Factor01)
fprintf(' SECTION 2 - Steel Section Ductility Ratio : %5.2f\n',Ductility_Rito02)
fprintf(' SECTION 2 - Steel Section Over Strength Factor : %5.2f\n',Over_Strength_Factor02)
fprintf(' SECTION 3 - Steel Section Ductility Ratio : %5.2f\n',Ductility_Rito03)
fprintf(' SECTION 3 - Steel Section Over Strength Factor : %5.2f\n',Over_Strength_Factor03)
fprintf(' SECTION 4 - Steel Section Ductility Ratio : %5.2f\n',Ductility_Rito04)
fprintf(' SECTION 4 - Steel Section Over Strength Factor : %5.2f\n',Over_Strength_Factor04)
fprintf('+-------------------------------------------------------+\n')
%%%  print time of computation
totaltime = cputime - starttime;
fprintf('\nTotal time (s): %7.4f \n\n',totaltime)
%% imaging
figure(1)
IMAGE1=imread('FiberSteelSectionMomentCurvature-image.jpg');
image(IMAGE1);axis image;axis off;
figure(2)
IMAGE1=imread('FiberSteelSectionMomentCurvatureMix.jpg');
image(IMAGE1);axis image;axis off;
%% Plot
figure(3)
p1=plot(I01,DX01,I02,DX02,'r--',I03,DX03,'g--',I04,DX04,'--');grid on;set(p1,'LineWidth',2);
xlabel('increment');ylabel('Residual');
legend('Section 1','Section 2','Section 3','Section 4','Location','NorthEastOutside');
title('Residual-increment diagram','color','b');
figure(4)
p1=plot(I01,IT01,I02,IT02,'r--',I03,IT03,'g--',I04,IT04,'--');grid on;set(p1,'LineWidth',3);
xlabel('increment');ylabel('Iteration');
legend('Section 1','Section 2','Section 3','Section 4','Location','NorthEastOutside');
title('Iteration-increment diagram','color','b');
figure(5)
plot(ES01,XX01,ES02,XX02,'r--',ES03,XX03,'g--',ES04,XX04,'--','LineWidth',3);
xlabel('Top steel strain (kN/mm^2)');ylabel('Neuteral axis (mm)')
title(' Ploting of the Top steel strain and Neuteral axis ','Color','b');grid on
legend('Section 1','Section 2','Section 3','Section 4','Location','NorthEastOutside');
figure(6)
plot(ES01,CUR01,ES02,CUR02,'r--',ES03,CUR03,'g--',ES04,CUR04,'--','LineWidth',3);
xlabel('Top concrete strain (mm/mm)');ylabel('Curvature (1/m)')
title(' Ploting of the Top concrete strain and Curvature ','Color','b');grid on
legend('Section 1','Section 2','Section 3','Section 4','Location','NorthEastOutside');
figure(7)
plot(CUR01,EI01,CUR02,EI02,'r--',CUR03,EI03,'g--',CUR04,EI04,'--','LineWidth',3)
title('#  Flextural Rigidity(EI)-CURVATURE DIAGRAM #','Color','b'); 
xlabel('CURVATURE (1/m)');ylabel('EI (kN.m^2)');grid on;
legend('Section 1','Section 2','Section 3','Section 4','Location','NorthEastOutside');
figure(8)
plot(XX01,EI01,XX02,EI02,'r--',XX03,EI03,'g--',XX04,EI04,'--','LineWidth',3)
title('# Neuteral axis-Flextural Rigidity(EI) diagram #','Color','b'); 
xlabel('Neuteral axis (mm)');ylabel('EI (kN.m^2)');grid on;
legend('Section 1','Section 2','Section 3','Section 4','Location','NorthEastOutside');
figure(9)
plot(Cur01,Mom01,Cur02,Mom02,'r--',Cur03,Mom03,'g--',Cur04,Mom04,'--','LineWidth',3)
title('# STEEL SECTION MOMENT CURVATURE DIAGRAM OF ALL SECTIONS #','Color','b'); 
xlabel('CURVATURE (1/m)');ylabel('MOMENT (kN-m)');grid on;
legend('Section 1','Section 2','Section 3','Section 4','Location','NorthEastOutside');
figure(10)
bar([Elastic_EI01,Elastic_EI02,Elastic_EI03,Elastic_EI04])
xlabel('Section Num(i)');ylabel('Elastic Flextural Rigidity EI_e_l_a (kN.m^2)');grid(gca,'minor');
figure(11)
bar([Plastic_EI01,Plastic_EI02,Plastic_EI03,Plastic_EI04])
xlabel('Section Num(i)');ylabel('Plastic Flextural Rigidity EI_p_l_a (kN.m^2)');grid(gca,'minor');
figure(12)
bar([Over_Strength_Factor01,Over_Strength_Factor02,Over_Strength_Factor03,Over_Strength_Factor04])
xlabel('Section Num(i)');ylabel('Over Strength Factor (Mu/My)(kN.m^2)');grid(gca,'minor');
%% Output Data to .txt file
fid = fopen('FiberSteelSectionMomentCurvatureMix-OutPut.txt','w');
fprintf(fid,'*************************************************************\n');
fprintf(fid,'*                  >> IN THE NAME OF GOD <<                 *\n');
fprintf(fid,'*      Moment-Curvature Analysis of 4 steel sections        *\n');
fprintf(fid,'*-----------------------------------------------------------*\n');
fprintf(fid,'*     This program is written by salar delavar ghashghaei   *\n');  
fprintf(fid,'*          E-mail:salar.d.ghashghaei@gmail.com              *\n');
fprintf(fid,'*-----------------------------------------------------------*\n');
fprintf(fid,'*Unit: Newton-Milimeter                                     *\n');
fprintf(fid,'*Given:Section Properties , Steel Section properties ,      *\n');
fprintf(fid,'*Calculate: Moment-Curavture                                *\n');
fprintf(fid,'*************************************************************\n');
fprintf(fid,'*        _     __                 __                        *\n');   
fprintf(fid,'*        ^    |  |               |  |                       *\n');
fprintf(fid,'*        |    |  |               |  |                       *\n');
fprintf(fid,'*             |  |               |  |                       *\n');
fprintf(fid,'*             |  |_______________|  |  +                    *\n');
fprintf(fid,'*       bf    |   _______________   |  tw                   *\n');
fprintf(fid,'*             |  |               |  |  +                    *\n');
fprintf(fid,'*        |    |  |               |  |                       *\n');
fprintf(fid,'*        v    |  |               |  |                       *\n');
fprintf(fid,'*        -     ---               ---                        *\n');
fprintf(fid,'*             |<-    2*tf+ hw     ->|                       *\n'); 
fprintf(fid,'*    X                                                      *\n');
fprintf(fid,'*    ^                                                      *\n');
fprintf(fid,'*    |             (Moment - Curvature along X axis)        *\n');
fprintf(fid,'*    |                                                      *\n');
fprintf(fid,'*    +----> Y                                               *\n');
fprintf(fid,'*************************************************************\n');
fprintf(fid,'\n');
fprintf(fid,'+==============================+\n');
fprintf(fid,'=          SECTION 1           =\n');
fprintf(fid,'     Curvature    Moment        \n');
fprintf(fid,'       (1/m)      (kN.m)        \n');
fprintf(fid,'--------------------------------\n');
fprintf(fid,'%10.5f %10.3f\n',[Cur01;Mom01]);
fprintf(fid,'+==============================+\n');
fprintf(fid,'=          SECTION 2           =\n');
fprintf(fid,'     Curvature    Moment        \n');
fprintf(fid,'       (1/m)      (kN.m)        \n');
fprintf(fid,'--------------------------------\n');
fprintf(fid,'%10.5f %10.3f\n',[Cur02;Mom02]);
fprintf(fid,'+==============================+\n');
fprintf(fid,'=          SECTION 3           =\n');
fprintf(fid,'     Curvature    Moment        \n');
fprintf(fid,'       (1/m)      (kN.m)        \n');
fprintf(fid,'--------------------------------\n');
fprintf(fid,'%10.5f %10.3f\n',[Cur03;Mom03]);
fprintf(fid,'+==============================+\n');
fprintf(fid,'=          SECTION 4           =\n');
fprintf(fid,'     Curvature    Moment        \n');
fprintf(fid,'       (1/m)      (kN.m)        \n');
fprintf(fid,'--------------------------------\n');
fprintf(fid,'%10.5f %10.3f\n',[Cur04;Mom04]);
fprintf(fid,'+======================================================================================================================+\n');
fprintf(fid,'=     Section 1 curve fitted  =    Section 2 curve fitted  =     Section 1 curve fitted  =     Section 1 curve fitted  =\n');
fprintf(fid,'     Curvature    Moment          Curvature    Moment            Curvature    Moment          Curvature    Moment       \n');
fprintf(fid,'       (1/m)      (kN.m)            (1/m)      (kN.m)              (1/m)      (kN.m)            (1/m)      (kN.m)       \n');
fprintf(fid,'------------------------------------------------------------------------------------------------------------------------\n');
fprintf(fid,'%10.4f %10.3f %20.4f %10.3f %20.4f %10.3f %20.4f %10.3f\n',[X01;Y01;X02;Y02;X03;Y03;X04;Y04]);
fprintf(fid,'+======================================================================================================================+\n');
fprintf(fid,'\n');
fprintf(fid,'+================================ SECTION 1 ===================================+\n');
fprintf(fid,' Increment   Top strain   Neuteral axis(x)    Curvature   Flextural Rigidity(EI)\n');
fprintf(fid,'================================================================================\n');
fprintf(fid,'  (i)           (1)           (mm)              (1/m)            (kN.m^2)       \n');
fprintf(fid,'--------------------------------------------------------------------------------\n');
fprintf(fid,'%4.0f %17.5f %13.2f %17.6f %17.2f\n',[I01;ES01;XX01;CUR01;EI01]);
fprintf(fid,'+================================ SECTION 2 ===================================+\n');
fprintf(fid,' Increment   Top strain   Neuteral axis(x)    Curvature   Flextural Rigidity(EI)\n');
fprintf(fid,'================================================================================\n');
fprintf(fid,'  (i)           (1)           (mm)              (1/m)            (kN.m^2)       \n');
fprintf(fid,'--------------------------------------------------------------------------------\n');
fprintf(fid,'%4.0f %17.5f %13.2f %17.6f %17.2f\n',[I02;ES02;XX02;CUR02;EI02]);
fprintf(fid,'+================================ SECTION 3 ===================================+\n');
fprintf(fid,' Increment   Top strain   Neuteral axis(x)    Curvature   Flextural Rigidity(EI)\n');
fprintf(fid,'================================================================================\n');
fprintf(fid,'  (i)           (1)           (mm)              (1/m)            (kN.m^2)       \n');
fprintf(fid,'--------------------------------------------------------------------------------\n');
fprintf(fid,'%4.0f %17.5f %13.2f %17.6f %17.2f\n',[I03;ES03;XX03;CUR03;EI03]);
fprintf(fid,'+================================ SECTION 4 ===================================+\n');
fprintf(fid,' Increment   Top strain   Neuteral axis(x)    Curvature   Flextural Rigidity(EI)\n');
fprintf(fid,'================================================================================\n');
fprintf(fid,'  (i)           (1)           (mm)              (1/m)            (kN.m^2)       \n');
fprintf(fid,'--------------------------------------------------------------------------------\n');
fprintf(fid,'%4.0f %17.5f %13.2f %17.6f %17.2f\n',[I04;ES04;XX04;CUR04;EI04]);
fprintf(fid,'+===============================================================================+\n');
fclose(fid);