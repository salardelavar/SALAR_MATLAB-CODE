%***********************************************************%
%                  >> IN THE NAME OF GOD <<                 %
%   Moment-Curvature Analysis of I steel section            %
%   with Plate on Flange                                    %
%-----------------------------------------------------------%
%     This program is written by salar delavar ghashghaei   %  
%          E-mail:salar.d.ghashghaei@gmail.com              %
%-----------------------------------------------------------%
%Unit: Newton-Milimeter                                     %
%Given: Steel section properties                            %
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
%    Y                                                      %
%    ^                                                      %
%    |             (Moment - Curvature along X axis)        %
%    |                                                      %
%    +----> X                                               %
%***********************************************************%
clear all;close all;clc;
%% Section Properties
tf1=9.2;% [mm] I section thickness on Top flange
bf1=110;% [mm] I section width on Top flange
tw=5.9;% [mm] I section thickness of Web
ptf1=10;% [mm] Plate section thickness on Top flange :. Note: ptf1 less than bf1
pbf1=70;% [mm] Plate section width on Top flange
hw=201.6;% [mm] Height of web
tf2=tf1;% [mm] I section thickness on Bottom flange
bf2=bf1;% [mm] I section width on Bottom flange
ptf2=ptf1;% [mm] Plate section thickness on Bottom flange
pbf2=pbf1;% [mm] Plate section width on Bottom flange
h=ptf1+ptf2+tf1+tf2+hw;% [mm] Height of Section
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
x=.25*h;% initial guess of Neuteral axis
%%% monitor cpu time
starttime = cputime;
%% ------------------ Newton Method Procedure - Section 1 -----------------------%
disp('#############');
disp('# SECTION 1 #');
disp('#############');
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
c01=[c1 c2 c3 c4 c5];
ES01=[.2*ey .4*ey .6*ey .8*ey ey .2*esh .4*esh .6*esh .8*esh esh .2*esu .4*esu .6*esu .8*esu esu];q01=size(ES01,2);
% ES01=ey:.001:esu;q01=size(ES01,2);
for j=1:q01;
    eS=ES01(j);
    it = 0; % initialize iteration count
    residual = 100; % initialize residual
    while (residual > tolerance)
        for z=1:N;% in this step: steel section force for each fiber is calculated
        es=eS*(x-c01(z))/x;            
        if and(z>=1,z<=roundn(N/10,0)) % Top Plate of flange
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
     Fs(z)=pbf1*R1*ptf1*fs;% force on top flange
     Fstan(z)=pbf1*R1*ptf1*fstan;% tangent steel force
        elseif and(z>roundn(N/10,0),z<=roundn(2*N/10,0)) % Top flange of I
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
     Fs(z)=bf1*R2*tf1*fs;% force on top flange
     Fstan(z)=bf1*R2*tf1*fstan;% tangent steel force
        elseif and(z>roundn(2*N/10,0),z<=roundn(8*N/10,0)) % Middle web
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
     Fs(z)=tw*R3*hw*fs;% force on web
     Fstan(z)=tw*R3*hw*fstan;% tangent steel force
     elseif and(z>roundn(8*N/10,0),z<=roundn(9*N/10,0)) % Below flange of I
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
      Fs(z)=bf2*R4*tf2*fs;% force on bottom flange
      Fstan(z)=bf2*R4*tf2*fstan;% tangent steel force
      elseif and(z>roundn(9*N/10,0),z<=N) % below Plate of flange
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
      Fs(z)=pbf2*R5*ptf2*fs;% force on bottom flange
      Fstan(z)=pbf2*R5*ptf2*fstan;% tangent steel force
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
if roundn(eS,-5) == esu;fprintf('\n      ## Strain Reached to Ultimate Strain: %1.4f ## \n\n',eS);end
Cur01=[0 Cur01];Mom01=[0 Mom01];
s01=size(Cur01,2);for i=1:s01-1;EI01(i)=(Mom01(i+1)-Mom01(i))/(Cur01(i+1)-Cur01(i));end % Flextural Rigidity
%% Steel Section bilinear fitting - Section 1
SIZE01=size(Mom01,2);
for i=1:SIZE01-1;
    hh01(i) = Cur01(i+1)-Cur01(i);
    Aa01(i)=(Mom01(i)+Mom01(i+1))*0.5*hh01(i);
end
Area01=sum(Aa01);k001 =Mom01(2)/Cur01(2);
fiy01 = (Mom01(i+1)*max(Cur01)*0.5-Area01)/(Mom01(i+1)*0.5 - k001*max(Cur01)*0.5);
My01 = k001*fiy01;
X01 = [0 fiy01 max(Cur01)];Y01 = [0 My01 Mom01(i+1)];
%% Elastic-Perfect Plastic Steel Section Moment-Curvature Analysis (simple)- Section 1
%  in this simple analysis : bf1=bf2 & tf1=tf2
Ie1=((tw*hw^3)/12)+2*(((bf1*tf1^3)/12)+(bf1*tf1)*(0.5*h-0.5*tf1-ptf1)^2)+2*(((pbf1*ptf1^3)/12)+(pbf1*ptf1)*(0.5*h-0.5*ptf1)^2);
phiY1=ey/(.5*h);phiU1=esu/(.5*h);
My1=(fy*Ie1)/(.5*h); % Yeild Moment
Mu1=fy*[2*(bf1*tf1*(.5*h-.5*tf1)+pbf1*ptf1*(.5*h-.5*ptf1))+(tw*hw^2)/4];  % Ultimate MomentEIe1=(My1/phiY1)*1e-9;
EIe1=(My1/phiY1)*1e-9;% Elastic Flextural Rigidity
EIp1=((Mu1-My1)/(phiU1-phiY1))*1e-9;% Plastic Flextural Rigidity
phi1=[0;phiY1;phiU1]*1000;Mi1=[0;My1;Mu1]*10^-6;
SD1=phiU1/phiY1;a1=Mu1/My1;
%% ------------------ Newton Method Procedure - Section 2 -----------------------%
disp('#############');
disp('# SECTION 2 #');
disp('#############');
clear c01 R1 R2 R3 R4 R5 Fs Fstan 
R1=(1/N);R2=(1/N);R3=(1/N);R4=(1/N);R5=(1/N);
for k=1:N;c01(k)=(.5*bf1-.5*pbf1)+(.5*R1+R1*(k-1))*pbf1;end % distance of each top Plate flange Fiber 
for k=1:N;c02(k)=(.5*R2+R2*(k-1))*bf1;end % distance of each top I flange Fiber
for k=1:N;c03(k)=(.5*bf1-.5*tw)+(.5*R3+R3*(k-1))*tw;end % distance of each web Fiber
for k=1:N;c04(k)=(.5*R4+R4*(k-1))*bf2;end % distance of each bottom I flange Fiber
for k=1:N;c05(k)=(.5*bf2-.5*pbf2)+(.5*R5+R5*(k-1))*pbf2;end % distance of each bottom Plate flange Fiber
ES02=[.2*ey .4*ey .6*ey .8*ey ey .2*esh .4*esh .6*esh .8*esh esh .2*esu .4*esu .6*esu .8*esu esu];q02=size(ES02,2);
% ES02=ey:.001:esu;q02=size(ES02,2);
x=.25*h;% initial guess of Neuteral axis
for j=1:q02;
    eS1=ES02(j)*.5*(bf1+pbf1)/bf1;
    eS2=ES02(j);
    eS3=ES02(j)*.5*(bf1+tw)/bf1;
    eS4=ES02(j);
    eS5=ES02(j)*.5*(bf2+pbf2)/bf1;
    it = 0; % initialize iteration count
    residual = 100; % initialize residual
    while (residual > tolerance)
        for z=1:N;% flange plate
        es=eS1*(x-c01(z))/x;            
    %---------------- Fs -------------%
       if and(es>0,es<=ey)
        fs=Es*es;
        fstan=(Es*eS1*c01(z))/(x)^2;
       elseif and(es<0,es>=(-ey))
        fs=Es*es;
        fstan=(Es*eS1*c01(z))/(x)^2;
       elseif  and(es>ey,es<=esh)
        fs=fy;
        fstan=0;
       elseif and(es<(-ey),es>=(-esh))
        fs=-fy;
        fstan=0;
        elseif  and(es>esh,es<=esu)
        fs=fu-(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS1*c01(z)*(fu-fy)*(((eS1*c01(z))/x)+esu-eS1))/(x^2*(esu-esh)^2);
       elseif and(es<(-esh),es>=(-esu))
        fs=-fu+(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS1*c01(z)*(fu-fy)*(((eS1*c01(z))/x)+esu-eS1))/(x^2*(esu-esh)^2);
       elseif or(es>esu,es<(-esu))
        fs=0;
        fstan=0;
       end
     Fs1(z)=pbf1*R1*ptf1*fs;% force on top flange
     Fstan1(z)=pbf1*R1*ptf1*fstan;% tangent steel force
     fs1(z)=fs;
        end
        for z=1:N;% flange of I
        es=eS2*(x-c02(z))/x;            
    %---------------- Fs -------------%
       if and(es>0,es<=ey)
        fs=Es*es;
        fstan=(Es*eS2*c02(z))/(x)^2;
       elseif and(es<0,es>=(-ey))
        fs=Es*es;
        fstan=(Es*eS2*c02(z))/(x)^2;
       elseif  and(es>ey,es<=esh)
        fs=fy;
        fstan=0;
       elseif and(es<(-ey),es>=(-esh))
        fs=-fy;
        fstan=0;
        elseif  and(es>esh,es<=esu)
        fs=fu-(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS2*c02(z)*(fu-fy)*(((eS2*c02(z))/x)+esu-eS2))/(x^2*(esu-esh)^2);
       elseif and(es<(-esh),es>=(-esu))
        fs=-fu+(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS2*c02(z)*(fu-fy)*(((eS2*c02(z))/x)+esu-eS2))/(x^2*(esu-esh)^2);
       elseif or(es>esu,es<(-esu))
        fs=0;
        fstan=0;
       end
     Fs2(z)=bf1*R2*tf1*fs;% force on top flange
     Fstan2(z)=bf1*R2*tf1*fstan;% tangent steel force
     fs2(z)=fs;
        end
        for z=1:N;% web of I
        es=eS3*(x-c03(z))/x;            
    %---------------- Fs -------------%
       if and(es>0,es<=ey)
        fs=Es*es;
        fstan=(Es*eS3*c03(z))/(x)^2;
       elseif and(es<0,es>=(-ey))
        fs=Es*es;
        fstan=(Es*eS3*c03(z))/(x)^2;
       elseif  and(es>ey,es<=esh)
        fs=fy;
        fstan=0;
       elseif and(es<(-ey),es>=(-esh))
        fs=-fy;
        fstan=0;
        elseif  and(es>esh,es<=esu)
        fs=fu-(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS3*c03(z)*(fu-fy)*(((eS3*c03(z))/x)+esu-eS3))/(x^2*(esu-esh)^2);
       elseif and(es<(-esh),es>=(-esu))
        fs=-fu+(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS3*c03(z)*(fu-fy)*(((eS3*c03(z))/x)+esu-eS3))/(x^2*(esu-esh)^2);
       elseif or(es>esu,es<(-esu))
        fs=0;
        fstan=0;
       end
     Fs3(z)=hw*R3*tw*fs;% force on top flange
     Fstan3(z)=hw*R3*tw*fstan;% tangent steel force
     fs3(z)=fs;
        end
        for z=1:N;% flange of I
        es=eS4*(x-c04(z))/x;            
    %---------------- Fs -------------%
       if and(es>0,es<=ey)
        fs=Es*es;
        fstan=(Es*eS4*c04(z))/(x)^2;
       elseif and(es<0,es>=(-ey))
        fs=Es*es;
        fstan=(Es*eS4*c04(z))/(x)^2;
       elseif  and(es>ey,es<=esh)
        fs=fy;
        fstan=0;
       elseif and(es<(-ey),es>=(-esh))
        fs=-fy;
        fstan=0;
        elseif  and(es>esh,es<=esu)
        fs=fu-(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS4*c04(z)*(fu-fy)*(((eS4*c04(z))/x)+esu-eS4))/(x^2*(esu-esh)^2);
       elseif and(es<(-esh),es>=(-esu))
        fs=-fu+(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS4*c04(z)*(fu-fy)*(((eS4*c04(z))/x)+esu-eS4))/(x^2*(esu-esh)^2);
       elseif or(es>esu,es<(-esu))
        fs=0;
        fstan=0;
       end
     Fs4(z)=bf2*R4*tf2*fs;% force on top flange
     Fstan4(z)=bf2*R4*tf2*fstan;% tangent steel force
     fs4(z)=fs;
        end
        for z=1:N;% flange plate
        es=eS5*(x-c05(z))/x;            
    %---------------- Fs -------------%
       if and(es>0,es<=ey)
        fs=Es*es;
        fstan=(Es*eS5*c05(z))/(x)^2;
       elseif and(es<0,es>=(-ey))
        fs=Es*es;
        fstan=(Es*eS5*c05(z))/(x)^2;
       elseif  and(es>ey,es<=esh)
        fs=fy;
        fstan=0;
       elseif and(es<(-ey),es>=(-esh))
        fs=-fy;
        fstan=0;
        elseif  and(es>esh,es<=esu)
        fs=fu-(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS5*c05(z)*(fu-fy)*(((eS5*c05(z))/x)+esu-eS5))/(x^2*(esu-esh)^2);
       elseif and(es<(-esh),es>=(-esu))
        fs=-fu+(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS5*c05(z)*(fu-fy)*(((eS5*c05(z))/x)+esu-eS5))/(x^2*(esu-esh)^2);
       elseif or(es>esu,es<(-esu))
        fs=0;
        fstan=0;
       end
     Fs5(z)=pbf2*R5*ptf2*fs;% force on top flange
     Fstan5(z)=pbf2*R5*ptf2*fstan;% tangent steel force
     fs5(z)=fs;
        end
    %----------------------------------%
    FsT1=sum(Fs1);FsT2=sum(Fs2);FsT3=sum(Fs3);FsT4=sum(Fs4);FsT5=sum(Fs5);
    A=FsT1+FsT2+FsT3+FsT4+FsT5;
    FsT1_tan=sum(Fstan1);FsT2_tan=sum(Fstan2);FsT3_tan=sum(Fstan3);FsT4_tan=sum(Fstan4);FsT5_tan=sum(Fstan5);
    A_tan=FsT1_tan+FsT2_tan+FsT3_tan+FsT4_tan+FsT5_tan;
    dx = A_tan^-1 *(-A);
    residual = abs(dx); % evaluate residual
        it = it + 1; % increment iteration count
        x = x+dx; % update x
        if it == itermax % stop the the analysis of this step please of Convergence
          fprintf('(-)Increment %1.0f : trail iteration reached to Ultimate %1.0f - strain: %1.6f - error: [%1.2f]\n',j,it,eS1,A)
             disp('    ## The solution for this step is not converged. Please check your model ##') 
            break
        end
    end
    if it == itermax;break;end % stop the analysis at all because last Convergence     
     for k=1:N;cc01(k)=x-c01(k);end %distance of each steel fiber from Neuteral axis of section 1
     for k=1:N;cc02(k)=x-c02(k);end %distance of each steel fiber from Neuteral axis of section 2
     for k=1:N;cc03(k)=x-c03(k);end %distance of each steel fiber from Neuteral axis of section 3
     for k=1:N;cc04(k)=x-c04(k);end %distance of each steel fiber from Neuteral axis of section 4
     for k=1:N;cc05(k)=x-c05(k);end %distance of each steel fiber from Neuteral axis of section 5

        if it < itermax% iteration control
        fprintf('(+)Increment %1.0f : It is converged in %1.0f iterations - strain: %1.5f - x: %1.2f - Phi: %1.5f - Moment: %1.2f\n',j,it,eS2,x,(eS2/x)*1000,(Fs1*cc01'+Fs2*cc02'+Fs3*cc03'+Fs4*cc04'+Fs5*cc05')*10^-6)
        end
 if and(eS2>0,eS2<ey);ESS02(j)=eS2;TUSS02(j)=Es*eS2;elseif and(eS2>=ey,eS2<esh);ESS02(j)=eS2;TUSS02(j)=fy;elseif and(eS2>=esh,eS2<esu);ESS02(j)=eS2;TUSS02(j)=fu-(fu-fy)*((esu-abs(eS2))/(esu-esh))^2;end% Top Steel section compression stress        
% Calculate Moment and Curavture
Cur02(j)=(eS2/x)*1000;XX02(j)=x;AA02(j)=A;
Mom02(j)=(Fs1*cc01'+Fs2*cc02'+Fs3*cc03'+Fs4*cc04'+Fs5*cc05')*10^-6;%1.0421
CUR02(j)=Cur02(j);I02(j)=j;IT02(j)=it;DX02(j)=dx;
end
if roundn(eS2,-5) == esu;fprintf('\n      ## Strain Reached to Ultimate Strain: %1.4f ## \n\n',eS2);end
Cur02=[0 Cur02];Mom02=[0 Mom02];
s02=size(Cur02,2);for i=1:s02-1;EI02(i)=(Mom02(i+1)-Mom02(i))/(Cur02(i+1)-Cur02(i));end % Flextural Rigidity
%% Steel Section bilinear fitting - Section 2
SIZE02=size(Mom02,2);
for i=1:SIZE02-1;
    hh02(i) = Cur02(i+1)-Cur02(i);
    Aa02(i)=(Mom02(i)+Mom02(i+1))*0.5*hh02(i);
end
Area02=sum(Aa02);k002 =Mom02(2)/Cur02(2);
fiy02 = (Mom02(i+1)*max(Cur02)*0.5-Area02)/(Mom02(i+1)*0.5 - k002*max(Cur02)*0.5);
My02 = k002*fiy02;
X02 = [0 fiy02 max(Cur02)];Y02 = [0 My02 Mom02(i+1)];
%% Elastic-Perfect Plastic Steel Section Moment-Curvature Analysis (simple)- Section 2
%  in this simple analysis : bf1=bf2 & tf1=tf2
Ie2=((hw*tw^3)/12)+(2*(tf1*bf1^3)/12)+(2*(ptf1*pbf1^3)/12);
phiY2=ey/(.5*bf1);phiU2=esu/(.5*bf1);
My2=(fy*Ie2)/(.5*bf1); % Yeild Moment
Mu2=fy*[(tw^2*hw)/4+2*(bf1^2*tf1)/4+2*(pbf1^2*ptf1)/4];  % Ultimate Moment
EIe2=(My2/phiY2)*1e-9;% Elastic Flextural Rigidity
EIp2=((Mu2-My2)/(phiU2-phiY2))*1e-9;% Plastic Flextural Rigidity
phi2=[0;phiY2;phiU2]*1000;Mi2=[0;My2;Mu2]*10^-6;
SD2=phiU2/phiY2;a2=Mu2/My2;
disp('+==========================================================+');
disp('=    Section 1 curve fitted   =   Section 2 curve fitted   =');
disp('     Curvature    Moment          Curvature    Moment       ');
disp('       (1/m)      (kN.m)            (1/m)      (kN.m)       ');
disp('------------------------------------------------------------');
disp([X01' Y01' X02' Y02']);
disp('+======================= Simple ===========================+');
disp('=         Section 1         =          Section 2           =');
disp('     Curvature    Moment          Curvature    Moment       ');
disp('       (1/m)      (kN.m)            (1/m)      (kN.m)       ');
disp('------------------------------------------------------------');
disp([phi1 Mi1 phi2 Mi2]);
disp('+==========================================================+');
%% EI and Ductility_Rito of Steel Section
Elastic_EI01=Y01(2)/X01(2);
Plastic_EI01=(Y01(3)-Y01(2))/(X01(3)-X01(2));
Ductility_Rito01=X01(3)/X01(2);
Over_Strength_Factor01=Y01(3)/Y01(2);
Elastic_EI02=Y02(2)/X02(2);
Plastic_EI02=(Y02(3)-Y02(2))/(X02(3)-X02(2));
Ductility_Rito02=X02(3)/X02(2);
Over_Strength_Factor02=Y02(3)/Y02(2);
Material_Ductility_Rito=esu/esh;
fprintf('+----------------------------------------------------------------+\n')
fprintf(' SECTION 1 - Elastic EI : %5.2f (kN.m^2)\n',Elastic_EI01)
fprintf(' SECTION 1 - Plastic EI : %5.2f (kN.m^2)\n',Plastic_EI01)
fprintf(' SECTION 1 - Elastic EI (Simple) : %5.2f (kN.m^2)\n',EIe1)
fprintf(' SECTION 1 - Plastic EI (Simple) : %5.2f (kN.m^2)\n',EIp1)
fprintf(' SECTION 2 - Elastic EI : %5.2f (kN.m^2)\n',Elastic_EI02)
fprintf(' SECTION 2 - Plastic EI : %5.2f (kN.m^2)\n',Plastic_EI02)
fprintf(' SECTION 2 - Elastic EI (Simple) : %5.2f (kN.m^2)\n',EIe2)
fprintf(' SECTION 2 - Plastic EI (Simple) : %5.2f (kN.m^2)\n',EIp2)
fprintf(' Steel Material Ductility Ratio : %5.2f\n',Material_Ductility_Rito)
fprintf(' SECTION 1 - Steel Section Ductility Ratio : %5.2f\n',Ductility_Rito01)
fprintf(' SECTION 1 - Steel Section Over Strength Factor : %5.2f\n',Over_Strength_Factor01)
fprintf(' SECTION 1 - Steel Section Ductility Ratio (Simple) : %5.2f\n',SD1)
fprintf(' SECTION 1 - Steel Section Over Strength Factor (Simple) : %5.2f\n',a1)
fprintf(' SECTION 2 - Steel Section Ductility Ratio : %5.2f\n',Ductility_Rito02)
fprintf(' SECTION 2 - Steel Section Over Strength Factor : %5.2f\n',Over_Strength_Factor02)
fprintf(' SECTION 2 - Steel Section Ductility Ratio (Simple) : %5.2f\n',SD2)
fprintf(' SECTION 2 - Steel Section Over Strength Factor (Simple) : %5.2f\n',a2)
fprintf('+----------------------------------------------------------------+\n')
%%%  print time of computation
totaltime = cputime - starttime;
fprintf('\nTotal time (s): %7.4f \n\n',totaltime)
%% imaging
figure(1)
IMAGE1=imread('FiberSteelSectionMomentCurvature-image.jpg');
image(IMAGE1);axis image;axis off;
figure(2)
IMAGE1=imread('FiberIPlateSteelSectionMomentCurvatureMix.jpg');
image(IMAGE1);axis image;axis off;
%% Plot
figure(3)
p1=plot(I01,DX01,I02,DX02,'r--');grid on;set(p1,'LineWidth',2);
xlabel('increment');ylabel('Residual');
title('Residual-increment diagram','color','b');
legend('Section 1','Section 2','Location','NorthEastOutside');grid on;
figure(4)
p1=plot(I01,IT01,I02,IT02,'r--');grid on;set(p1,'LineWidth',3);
xlabel('increment');ylabel('Iteration');
title('Iteration-increment diagram','color','b');
legend('Section 1','Section 2','Location','NorthEastOutside');grid on;
figure(5)
plot(ES01,XX01,ES02,XX02,'r--','LineWidth',3);xlabel('Top steel strain ');ylabel('Neuteral axis (mm)')
title(' Plotting of the Top steel strain and Neuteral axis ','Color','b');
legend('Section 1','Section 2','Location','NorthEastOutside');grid on;
figure(6)
plot(ES01,CUR01,ES02,CUR02,'r--','LineWidth',2);
xlabel('Top steel strain ');ylabel('Curvature (1/m)')
title(' Plotting of the Top steel strain and Curvature ','Color','b');
legend('Section 1','Section 2','Location','NorthEastOutside');grid on;
figure(7)
plot(CUR01,EI01,CUR02,EI02,'r--','LineWidth',3)
title('#  Flextural Rigidity(EI)-CURVATURE DIAGRAM #','Color','b'); 
xlabel('CURVATURE (1/m)');ylabel('EI (kN.m^2)');
legend('Section 1','Section 2','Location','NorthEastOutside');grid on;
figure(8)
plot(XX01,EI01,XX02,EI02,'r--','LineWidth',3)
title('# Neuteral Axis-Flextural Rigidity(EI) diagram #','Color','b'); 
xlabel('Neuteral axis (mm)');ylabel('EI (kN.m^2)');
legend('Section 1','Section 2','Location','NorthEastOutside');grid on;
figure(9)
plot(Cur01,Mom01,X01,Y01,'--',phi1,Mi1,'--',Cur02,Mom02,X02,Y02,'--',phi2,Mi2,'--','LineWidth',3)
title('# STEEL SECTION MOMENT CURVATURE DIAGRAM #','Color','b'); 
xlabel('CURVATURE (1/m)');ylabel('MOMENT (kN-m)')
legend('Section 1','Section 1 - bilinear fitted','Section 1 - Elastic-Perfect Plastic','Section 2','Section 2 - bilinear fitted','Section 2 - Elastic-Perfect Plastic','Location','NorthEastOutside');grid on;
%% Output Data to .txt file
fid = fopen('FiberIPlateSteelSectionMomentCurvatureMix-OutPut.txt','w');
fprintf(fid,'*************************************************************\n');
fprintf(fid,'*                  >> IN THE NAME OF GOD <<                 *\n');
fprintf(fid,'*   Moment-Curvature Analysis of I steel section            *\n');
fprintf(fid,'*   with Plate on Flange                                    *\n');
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
fprintf(fid,'*    Y                                                      *\n');
fprintf(fid,'*    ^                                                      *\n');
fprintf(fid,'*    |             (Moment - Curvature along X axis)        *\n');
fprintf(fid,'*    |                                                      *\n');
fprintf(fid,'*    +----> X                                               *\n');
fprintf(fid,'*************************************************************\n');
fprintf(fid,'\n');
fprintf(fid,'+==============================+\n');
fprintf(fid,'=            Sectin 1          =\n');
fprintf(fid,'     Curvature    Moment        \n');
fprintf(fid,'       (1/m)      (kN.m)        \n');
fprintf(fid,'--------------------------------\n');
fprintf(fid,'%10.5f %10.3f\n',[Cur01;Mom01]);
fprintf(fid,'+==============================+\n');
fprintf(fid,'=            Sectin 2          =\n');
fprintf(fid,'     Curvature    Moment        \n');
fprintf(fid,'       (1/m)      (kN.m)        \n');
fprintf(fid,'--------------------------------\n');
fprintf(fid,'%10.5f %10.3f\n',[Cur02;Mom02]);
fprintf(fid,'+==========================================================+\n');
fprintf(fid,'= Steel Section curve fitted  = Steel Section Simple Anal. =\n');
fprintf(fid,'     Curvature    Moment          Curvature    Moment       \n');
fprintf(fid,'       (1/m)      (kN.m)            (1/m)      (kN.m)       \n');
fprintf(fid,'------------------------------------------------------------\n');
fprintf(fid,'%10.4f %10.3f %20.4f %10.3f\n',[X01;Y01]);
fprintf(fid,'+==========================================================+\n');
fprintf(fid,'\n');
fprintf(fid,'+================================== Section 1 =================================+\n');
fprintf(fid,' Increment   Top strain   Neuteral axis(x)    Curvature   Flextural Rigidity(EI)\n');
fprintf(fid,'================================================================================\n');
fprintf(fid,'  (i)           (1)           (mm)              (1/m)            (kN.m^2)       \n');
fprintf(fid,'--------------------------------------------------------------------------------\n');
fprintf(fid,'%4.0f %17.5f %13.2f %17.6f %17.2f\n',[I01;ES01;XX01;CUR01;EI01]);
fprintf(fid,'+================================== Section 2 =================================+\n');
fprintf(fid,' Increment   Top strain   Neuteral axis(x)    Curvature   Flextural Rigidity(EI)\n');
fprintf(fid,'================================================================================\n');
fprintf(fid,'  (i)           (1)           (mm)              (1/m)            (kN.m^2)       \n');
fprintf(fid,'--------------------------------------------------------------------------------\n');
fprintf(fid,'%4.0f %17.5f %13.2f %17.6f %17.2f\n',[I02;ES02;XX02;CUR02;EI02]);
fprintf(fid,'================================================================================\n');
fclose(fid);