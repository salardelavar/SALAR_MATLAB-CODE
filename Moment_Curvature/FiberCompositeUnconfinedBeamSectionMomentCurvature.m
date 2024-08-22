%***********************************************************%
%                  >> IN THE NAME OF GOD <<                 %
%  Moment-Curvature Analysis of Composite Unconfined Beam   %
%  Section                                                  %
%-----------------------------------------------------------%
%     This program is written by salar delavar ghashghaei   %  
%          E-mail:salar.d.ghashghaei@gmail.com              %
%-----------------------------------------------------------%
%Unit: Newton-Milimeter                                     %
%Given:Section Properties , SteelSection properties ,       %
%Calculate: Moment-Curavture                                %
% Newton-Raphson Method : Tangent procedure                 %
%***********************************************************%
%   _    ______________________________________             %
%   |   |                                      |            %
%       |     #     #     #     #    #    #    |            %
%       |     #                           #    |            %
% Beff  |    As1   As2   As3   As4  As5  As6   |            %
%       |     #                           #    |            %
%   |   |     #     #     #     #    #    #    |            %
%   _   |______________________________________|            %
%       |<-               Teff               ->|            %
%       |<-d1->|                                            %
%       |<-  d2   ->|                                       %
%       |<-     d3      ->|                                 %
%       |<-        d4          ->|                          %
%       |<-            d5          ->|                      %
%       |<-               d6             >|                 %
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
%    X                                                      %
%    ^                                                      %
%    |             (Moment - Curvature along X axis)        %
%    |                                                      %
%    +----> Y                                               %
%***********************************************************%
clear all;close all;clc;
%% Section Properties
Beff=1000; % Width of deck
Teff=100; % Depth of deck
%As:  As1      As2
As=[785.375 785.375];
%d:d1  d2 
d=[30 70];
tf1=9.2;% [mm] I section thickness on Top flange
bf1=110;% [mm] I section width on Top flange
tw=5.9;% [mm] I section thickness of Web
hw=201.6;% [mm] Height of web
tf2=9.2;% [mm] I section thickness on Bottom flange
bf2=110;% [mm] I section width on Bottom flange
ptf2=10;% [mm] Plate section thickness on Bottom flange
pbf2=80;% [mm] Plate section width on Bottom flange
h=Teff+tf1+tf2+hw+ptf2;% [mm] Height of Section
Area = tf1*bf1 +Beff*Teff +tw*hw + tf2*bf2 +ptf2*pbf2;
%% Steel Section Properties
fy =240;% [N/mm^2] Yield strength of steel section
Es =2e5;% [N/mm^2] Modulus of elasticity of steel section
fu=1.5*fy;% Ultimate steel stress
ey=fy/Es;% Yield steel strain
esh=0.025;% Strain at steel strain-hardening
esu=0.35;% Ultimate steel strain
Esh=(fu-fy)/(esu-esh);
%% Concrete Properties
fc =25;% [N/mm^2] Unconfined concrete strength
ecu=0.004;% Ultimate concrete strain
Ec=5000*sqrt(fc);
ec0=(2*fc)/Ec;
fct=-0.7*sqrt(fc);% Concrete tension stress
ect1=(2*fct)/Ec;ect2=(2.625*fct)/Ec;ect3=(9.292*fct)/Ec;% Concrete tension strain
%% Steel Reinforcing Properties
fys =400;% Yield strength of steel reinforcing (N/mm^2)
Ess =2e5;% Modulus of elasticity of steel (N/mm^2)
fus=1.5*fys;% Ultimate steel stress
eys=fys/Ess;% Yeild steel strain
eshs=0.01;% Strain at steel strain-hardening
esus=0.09;% Ultimate steel strain
Eshs=(fus-fys)/(esus-eshs);
N=h*2;% Number of steel section Fiber 
itermax = 4000;% maximum number of iterations
tolerance = 1e-12; % specified tolerance for convergence
x=.25*h;% initial guess of Neuteral axis
%%% monitor cpu time
starttime = cputime;
%% ------------------ Newton Method Procedure ------------------------%
N1=2*N/10;% Number of concrete section  
N2=N/10;% Number of steel section top flange Fiber  
N3=5*N/10;% Number of steel section  web Fiber
N4=N/10;% Number of steel section bottom flange Fiber
N5=N/10;% Number of steel section bottom Plate flange Fiber 
An=size(As,2);R1=(1/N1);R2=(1/N2);R3=(1/N3);R4=(1/N4);R5=(1/N5);
for k=1:N1;c1(k)=(.5*R1+R1*(k-1))*Teff;end 
for k=1:N2;c2(k)=Teff+(.5*R2+R2*(k-1))*tf1;end % distance of each top flange Fiber 
for k=1:N3;c3(k)=Teff+tf1+(.5*R3+R3*(k-1))*hw;end % distance of each web Fiber
for k=1:N4;c4(k)=Teff+tf1+hw+(.5*R4+R4*(k-1))*tf2;end % distance of each bottom flange Fiber
for k=1:N5;c5(k)=Teff+tf1+hw+tf2+(.5*R5+R5*(k-1))*ptf2;end % distance of each bottom flange Fiber
c=[c1 c2 c3 c4 c5];
ES=[.1*abs(ect1) .33*abs(ect1) .67*abs(ect1) abs(ect1) .8*abs(ect2) .9*abs(ect2) abs(ect2) ...
    .33*abs(ect3) .5*abs(ect3) .6*abs(ect3) .67*abs(ect3) .7*abs(ect3) .75*abs(ect3) .8*abs(ect3)...
    .9*abs(ect3) abs(ect3) .35*ecu .4*ecu .5*ecu .6*ecu .7*ecu .8*ecu .9*ecu ecu];q=size(ES,2);
% ES=.5*abs(ect1):.0001:ecu;q=size(ES,2);
for j=1:q;
    eS=ES(j);
    it = 0; % initialize iteration count
    residual = 100; % initialize residual
    while (residual > tolerance)
        %---------------- As -------------%
        for u=1:An;
        ess=eS*(x-d(u))/x;
       if and(ess>=0,ess<eys)
        fsa(u)=Ess*ess;
        fstana(u)=(Ess*eS*d(u))/(x)^2;
       elseif and(ess<0,ess>(-eys))
        fsa(u)=Ess*ess;
        fstana(u)=(Ess*eS*d(u))/(x)^2;
       elseif  and(ess>=eys,ess<eshs)
        fsa(u)=fys;
        fstana(u)=0;
       elseif and(ess<=(-eys),ess>(-eshs))
        fsa(u)=-fys;
        fstana(u)=0;
        elseif  and(ess>=eshs,ess<esus)
        fsa(u)=fus-(fus-fys)*((esus-abs(ess))/(esus-eshs))^2;
        fstana(u)=(2*eS*d(u)*(fus-fys)*(((eS*d(u))/x)+esus-eS))/(x^2*(esus-eshs)^2);
       elseif and(ess<=(-eshs),ess>(-esus))
        fsa(u)=-fus+(fus-fys)*((esus-abs(ess))/(esus-eshs))^2;
        fstana(u)=(2*eS*d(u)*(fus-fys)*(((eS*d(u))/x)+esus-eS))/(x^2*(esus-eshs)^2);
       elseif or(ess>=esus,ess<=(-esus))
        fsa(u)=0;
        fstana(u)=0;
       end
        FSA(j,u)=fsa(u);
        FsA(u)=As(u)*fsa(u);
        FstanA(u)=As(u)*fstana(u);% tangent steel force
        if and(ess>=0,ess<ec0)% in this step: Unconfined concrete force in rebar area is omitted (F=As*fc)
        Cs=fc*((2*ess/ec0)-(ess/ec0)^2);
        Ctans=((2*fc)/(ec0^2*x^3))*(ec0*eS*d(u)*x-d(u)*eS^2*x+2*eS^2*d(u)^2);
      elseif and(ess>=ec0,ess<ecu)
        Cs=fc*(1-(0.15*(ess-ec0)/(ecu-ec0)));
        Ctans=-(3*eS*d(u)*fc)/(20*(ecu-ec0)*x^2);
      elseif ess >= ecu
         Cs=0; 
         Ctans=0;
      elseif and(ess<0,ess>=ect1)
        Cs=-0.5*Ec*ess;
        Ctans=-(0.5*Ec*d(u)*eS)/x^2;
      elseif and(ess<ect1,ess>=ect2)
        Cs=fct+(0.5*fct/(ect2-ect1))*(ess-ect1);
        Ctans=+(0.5*fct*eS*d(u))/((ect2-ect1)*x^2);
      elseif and(ess<ect2,ess>=ect3)
        Cs=-(.5*fct-(0.5*fct/(ect3-ect2))*(ess-ect2));
        Ctans=+(0.5*fct*eS*d(u))/((ect3-ect2)*x^2);
      elseif ess<ect3
        Cs=0;
        Ctans=0;
        end
        CS(u)=Cs;CtanS(u)=Ctans;
        Fss(u)=-As(u)*CS(u);
        Fstans(u)=-As(u)*CtanS(u);% tangent Minus of concrete force
        end
        for z=1:N;% in this step: steel section force for each fiber is calculated
        es=eS*(x-c(z))/x;            
        if and(z>=1,z<=roundn(2*N/10,0)) % TOP CONCRETE SECTION
        if and(es>=0,es<ec0)
        fs=fc*((2*es/ec0)-(es/ec0)^2);
        fstan=((2*fc)/(ec0^2*x^3))*(ec0*eS*c(z)*x-c(z)*eS^2*x+2*eS^2*c(z)^2);
      elseif and(es>=ec0,es<=ecu)
        fs=fc*(1-(0.15*(es-ec0)/(ecu-ec0)));
        fstan=-(3*eS*c(z)*fc)/(20*(ecu-ec0)*x^2);
      elseif es > ecu
         fs=0; 
         fstan=0;
      elseif and(es<0,es>=ect1)
        fs=0.5*Ec*es;
        fstan=(0.5*Ec*c(z)*eS)/x^2;
      elseif and(es<ect1,es>=ect2)
        fs=fct-(0.5*fct/(ect2-ect1))*(es-ect1);
        fstan=-(0.5*fct*eS*c(z))/((ect2-ect1)*x^2);
      elseif and(es<ect2,es>=ect3)
        fs=.5*fct-(0.5*fct/(ect3-ect2))*(es-ect2);
        fstan=-(0.5*fct*eS*c(z))/((ect3-ect2)*x^2);
      elseif es<ect3
        fs=0;
        fstan=0;
        end
      Fs(z)=Beff*R1*Teff*fs;
      Fstan(z)=Beff*R1*Teff*fstan;% tangent concrete force  
        elseif and(z>=roundn(2*N/10,0),z<=roundn(3*N/10,0)) % TOP PLATE FLANGE    
        %---------------- Fs -------------%
       if and(es>=0,es<=ey)
        fs=Es*es;
        fstan=(Es*eS*c(z))/(x)^2;
       elseif and(es<0,es>=(-ey))
        fs=Es*es;
        fstan=(Es*eS*c(z))/(x)^2;
       elseif  and(es>ey,es<=esh)
        fs=fy;
        fstan=0;
       elseif and(es<(-ey),es>=(-esh))
        fs=-fy;
        fstan=0;
        elseif  and(es>esh,es<=esu)
        fs=fu-(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c(z)*(fu-fy)*(((eS*c(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif and(es<(-esh),es>=(-esu))
        fs=-fu+(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c(z)*(fu-fy)*(((eS*c(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif or(es>esu,es<(-esu))
        fs=0;
        fstan=0;
       end
     Fs(z)=bf1*R2*tf1*fs;% force on top flange
     Fstan(z)=bf1*R2*tf1*fstan;% tangent steel force
        elseif and(z>roundn(3*N/10,0),z<=roundn(8*N/10,0)) % MIDDLE WEB
       if and(es>=0,es<=ey)
        fs=Es*es;
        fstan=(Es*eS*c(z))/(x)^2;
       elseif and(es<0,es>=(-ey))
        fs=Es*es;
        fstan=(Es*eS*c(z))/(x)^2;
       elseif  and(es>ey,es<=esh)
        fs=fy;
        fstan=0;
       elseif and(es<(-ey),es>=(-esh))
        fs=-fy;
        fstan=0;
        elseif  and(es>esh,es<=esu)
        fs=fu-(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c(z)*(fu-fy)*(((eS*c(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif and(es<(-esh),es>=(-esu))
        fs=-fu+(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c(z)*(fu-fy)*(((eS*c(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif or(es>esu,es<(-esu))
        fs=0;
        fstan=0;
       end
     Fs(z)=tw*R3*hw*fs;% force on web
     Fstan(z)=tw*R3*hw*fstan;% tangent steel force
        elseif and(z>roundn(8*N/10,0),z<=roundn(9*N/10,0)) % BOTTOM I FLANGE
       if and(es>=0,es<=ey)
        fs=Es*es;
        fstan=(Es*eS*c(z))/(x)^2;
       elseif and(es<0,es>=(-ey))
        fs=Es*es;
        fstan=(Es*eS*c(z))/(x)^2;
       elseif  and(es>ey,es<=esh)
        fs=fy;
        fstan=0;
       elseif and(es<(-ey),es>=(-esh))
        fs=-fy;
        fstan=0;
        elseif  and(es>esh,es<=esu)
        fs=fu-(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c(z)*(fu-fy)*(((eS*c(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif and(es<(-esh),es>=(-esu))
        fs=-fu+(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c(z)*(fu-fy)*(((eS*c(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif or(es>esu,es<(-esu))
        fs=0;
        fstan=0;
       end
      Fs(z)=bf2*R4*tf2*fs;% force on bottom flange
      Fstan(z)=bf2*R4*tf2*fstan;% tangent steel force
     elseif and(z>roundn(9*N/10,0),z<=N) % BOTTOM PLATE FLANGE
       if and(es>=0,es<=ey)
        fs=Es*es;
        fstan=(Es*eS*c(z))/(x)^2;
       elseif and(es<0,es>=(-ey))
        fs=Es*es;
        fstan=(Es*eS*c(z))/(x)^2;
       elseif  and(es>ey,es<=esh)
        fs=fy;
        fstan=0;
       elseif and(es<(-ey),es>=(-esh))
        fs=-fy;
        fstan=0;
        elseif  and(es>esh,es<=esu)
        fs=fu-(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c(z)*(fu-fy)*(((eS*c(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif and(es<(-esh),es>=(-esu))
        fs=-fu+(fu-fy)*((esu-abs(es))/(esu-esh))^2;
        fstan=(2*eS*c(z)*(fu-fy)*(((eS*c(z))/x)+esu-eS))/(x^2*(esu-esh)^2);
       elseif or(es>esu,es<(-esu))
        fs=0;
        fstan=0;
       end
      Fs(z)=pbf2*R5*ptf2*fs;% force on plate bottom flange
      Fstan(z)=pbf2*R5*ptf2*fstan;% tangent steel force
        end 
     Ss(z)=es;SS(j,z)=Ss(z);% sction Fiber Strain
     CFS(j,z)=fs;% sction Fiber Stress
        end
    %----------------------------------%
    FsTOTAL=sum(Fs);FsATOTAL=sum(FsA);FssTOTAL=sum(Fss);A=FsTOTAL+FsATOTAL+FssTOTAL;
    FsTOTAL_tan=sum(Fstan);FsATOTAL_tan=sum(FstanA);FssTOTAL_tan=sum(Fstans);A_tan=FsTOTAL_tan+FsATOTAL_tan+FssTOTAL_tan;
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
    for u=1:An;e(u)=x-d(u);end % distance of each rebar from neutral axis        
    if it == itermax;break;end % stop the analysis at all because last Convergence     
     for k=1:N;cc(k)=x-c(k);end %distance of each steel fiber from Neuteral axis
        if it < itermax% iteration control
        fprintf('(+)Increment %1.0f : It is converged in %1.0f iterations - strain: %1.5f - x: %1.2f - Phi: %1.5f - Moment: %1.2f\n',j,it,eS,x,(eS/x)*1000,(Fs*cc'+FsA*e'+Fss*e')*10^-6)
        end                        
TSCSs(j)=eS*(x-d(1))/x;TSCS(j)=fsa(1);% Top Steel compression strain-stress
BSCSs(j)=eS*(x-d(2))/x;BSCS(j)=fsa(2);% Bottom Steel compression strain-stress       
if and(eS>0,eS<ec0);ECC(j)=eS;TUCS(j)=fc*((2*eS/ec0)-(eS/ec0)^2);elseif and(eS>=ec0,eS<=ecu);ECC(j)=eS;TUCS(j)=fc*(1-(0.15*(eS-ec0)/(ecu-ec0)));end% Top Unconfined compression stress 
eST=eS*(x-Teff)/x;if and(eST>0,eST<ey);EST(j)=eST;TUST(j)=Es*eST;elseif and(eST>=ey,eST<esh);EST(j)=eST;TUST(j)=fy;elseif and(eST>=esh,eST<esu);EST(j)=eST;TUST(j)=fu-(fu-fy)*((esu-abs(eST))/(esu-esh))^2;elseif and(eST<0,eST>(-ey));EST(j)=eST;TUST(j)=Es*eST;elseif and(eST<=(-ey),eST>(-esh));EST(j)=eST;TUST(j)=-fy;elseif and(eST<=(-esh),eST>(-esu));EST(j)=eST;TUST(j)=-fu+(fu-fy)*((esu-abs(eST))/(esu-esh))^2;end% Top Steel section compression stress
eSB=eS*(x-h)/x;if and(eSB>0,eSB<ey);ESB(j)=eSB;TUSB(j)=Es*eSB;elseif and(eSB>=ey,eSB<esh);ESB(j)=eSB;TUSB(j)=fy;elseif and(eSB>=esh,eSB<esu);ESB(j)=eSB;TUSB(j)=fu-(fu-fy)*((esu-abs(eSB))/(esu-esh))^2;elseif and(eSB<0,eSB>(-ey));ESB(j)=eSB;TUSB(j)=Es*eSB;elseif and(eSB<=(-ey),eSB>(-esh));ESB(j)=eSB;TUSB(j)=-fy;elseif and(eSB<=(-esh),eSB>(-esu));ESB(j)=eSB;TUSB(j)=-fu+(fu-fy)*((esu-abs(eSB))/(esu-esh))^2;end% Bottom Steel section compression stress
% Calculate Moment and Curavture
Cur(j)=(eS/x)*1000;XX(j)=x;AA(j)=A;
Mom(j)=(Fs*cc'+FsA*e'+Fss*e')*10^-6;
CUR(j)=Cur(j);I(j)=j;IT(j)=it;DX(j)=dx;
end
Cur=[0 Cur];Mom=[0 Mom];
%s=size(Cur,2);for i=1:s-1;EI(i)=(Mom(i+1)-Mom(i))/(Cur(i+1)-Cur(i));end % Flextural Rigidity
s=size(Cur,2);for i=1:s-1;EI(i)=Mom(i)/Cur(i);end % Flextural Rigidity
if eS == ecu;fprintf('\n      ##  Concrete Strain Reached to Ultimate Strain: %1.4f ## \n\n',eS);end
%% Steel Section bilinear fitting
SIZE=size(Mom,2);
for i=1:SIZE-1;
    hh(i) = Cur(i+1)-Cur(i);
    Aa(i)=(Mom(i)+Mom(i+1))*0.5*hh(i);
end
Area=sum(Aa);k0 =Mom(2)/Cur(2);
fiy = (Mom(i+1)*max(Cur)*0.5-Area)/(Mom(i+1)*0.5 - k0*max(Cur)*0.5);
My = k0*fiy;
X = [0 fiy max(Cur)];Y = [0 My Mom(i+1)];
disp('+=============================+');
disp('=    Section curve fitted     =');
disp('     Curvature    Moment       ');
disp('       (1/m)      (kN.m)       ');
disp('-------------------------------');
disp([X' Y']);
disp('+==============================+');
%% EI and Ductility_Rito of Steel Section
Elastic_EI=Y(2)/X(2);
Plastic_EI=(Y(3)-Y(2))/(X(3)-X(2));
Ductility_Rito=X(3)/X(2);
Over_Strength_Factor=Y(3)/Y(2);
Material_Ductility_Rito=esu/esh;
fprintf('+------------------------------------------+\n')
fprintf(' Elastic EI : %5.2f (kN.m^2)\n',Elastic_EI)
fprintf(' Plastic EI : %5.2f (kN.m^2)\n',Plastic_EI)
fprintf(' Steel Material Ductility Ratio : %5.2f\n',Material_Ductility_Rito)
fprintf(' Steel Section Ductility Ratio : %5.2f\n',Ductility_Rito)
fprintf(' Steel Section Over Strength Factor : %5.2f\n',Over_Strength_Factor)
fprintf('+------------------------------------------+\n')
%%%  print time of computation
totaltime = cputime - starttime;
fprintf('\nTotal time (s): %7.4f \n\n',totaltime)
%% Moment-Curvature of B100x10-20fi10
CurP=[0
    5.455*10^-3
    6.234*10^-3
    1.065*10^-2
    1.714*10^-2
    2.338*10^-2
    2.935*10^-2
    4.182*10^-2
    5.61*10^-2
    7.818*10^-2];
MomP=[0
    154.33
    161.36
    179.66
    193.23
    200.37
    205.54
    207.28
    208.36
    209.37];
%% C++ outout FiberRecSectionMomentCurvatureLayerStrainComposite
Curvature=[0
1.74017e-006
3.46374e-006
5.15729e-006
7.57372e-006
1.02629e-005
1.30965e-005
1.61783e-005
1.93503e-005
2.25646e-005
2.57885e-005
2.92791e-005
3.32285e-005
3.72229e-005
4.11346e-005
4.49533e-005
4.87029e-005
5.24009e-005
5.6029e-005
5.95988e-005
6.31289e-005
6.66226e-005
7.00666e-005
7.34833e-005
7.68698e-005
8.02253e-005
];

Moment=[0
5.18369e+007
1.02761e+008
1.51913e+008
1.70448e+008
1.80204e+008
1.87357e+008
1.9277e+008
1.97385e+008
2.01467e+008
2.05133e+008
2.07621e+008
2.08646e+008
2.09378e+008
2.10016e+008
2.10599e+008
2.1112e+008
2.11583e+008
2.12024e+008
2.12437e+008
2.12817e+008
2.1317e+008
2.13509e+008
2.13827e+008
2.14126e+008
2.14411e+008
];
%% imaging
figure(1)
IMAGE1=imread('FiberCompositeUnconfinedBeamSectionMomentCurvature-image1.jpg');
image(IMAGE1);axis image;axis off;
figure(2)
IMAGE1=imread('FiberCompositeUnconfinedBeamSectionMomentCurvature-image2.jpg');
image(IMAGE1);axis image;axis off;
%% Plot
figure(3)
plot(TSCSs,TSCS,BSCSs,BSCS,'g--','LineWidth',3);xlabel('Steel strain (\epsilon_s)');ylabel('Steel stress (N/mm^2)')
title(' Plotting of the Top  and Bottom steel rebar strain-stress during the analysis','Color','b');
legend('Top steel rebar','Bottom steel rebar','Location','NorthEastOutside');grid on;
figure(4)
P1=plot(EST,TUST,ESB,TUSB,'r--');set(P1,'LineWidth',3);
xlabel('Strain');ylabel('Stress (kN/mm^2)')
legend('Top steel fiber','Bottom steel fiber','Location','NorthEastOutside');grid on;
title('Top and Bottom Steel Fiber Section Stress-Strain Diagram','color','b')
figure(5)
P1=plot(ECC,TUCS);set(P1,'LineWidth',3);
xlabel('Strain');ylabel('Stress (kN/mm^2)');grid on;
title('Top Concrete Fiber Section Stress-Strain Diagram','color','b')
figure(6)
p1=plot(I,DX,'b--');grid on;set(p1,'LineWidth',2);
xlabel('increment');ylabel('Residual');
title('Residual-increment diagram','color','b');
figure(7)
p1=plot(I,IT,'b--');grid on;set(p1,'LineWidth',3);
xlabel('increment');ylabel('Iteration');
title('Iteration-increment diagram','color','b');
figure(8)
plot(ES,XX,'--','LineWidth',3);xlabel('Top Concrete strain ');ylabel('Neuteral axis (mm)')
title(' Plotting of the Top concrete strain and Neuteral axis ','Color','b');grid on
figure(9)
plot(ES,CUR,'--','LineWidth',2);xlabel('Top concrete strain ');ylabel('Curvature (1/m)')
title(' Plotting of the Top concrete strain and Curvature ','Color','b');grid on
figure(10)
plot(SS(1,:),(h-c),SS(roundn(.25*q,0),:),(h-c),SS(roundn(.5*q,0),:),(h-c),SS(roundn(.75*q,0),:),(h-c),SS(roundn(q,0),:),(h-c),'g--','LineWidth',3);xlabel('Strain (mm/mm)');ylabel('Height (mm)')
title(' Plotting of the steel and concrete fiber strain and Heigh of section - Negative value: Tension - Positive value: Compression','Color','b');grid on
legend(['\phi= ',num2str(Cur(2))],['\phi= ',num2str(Cur(roundn(.25*q,0)))],['\phi= ',num2str(Cur(roundn(.5*q,0)))],['\phi= ',num2str(Cur(roundn(.75*q,0)))],['\phi= ',num2str(Cur(roundn(q+1,0)))],'Location','NorthEastOutside');
figure(11)
plot(CFS(1,:),(h-c),CFS(roundn(.25*q,0),:),(h-c),CFS(roundn(.5*q,0),:),(h-c),CFS(roundn(.75*q,0),:),(h-c),CFS(roundn(q,0),:),(h-c),'g--','LineWidth',3);xlabel('Stress (Mpa)');ylabel('Height (mm)')
title(' Plotting of the steel and concrete fiber stress and Heigh of section - Negative value: Tension - Positive value: Compression','Color','b');grid on
legend(['\phi= ',num2str(Cur(2))],['\phi= ',num2str(Cur(roundn(.25*q,0)))],['\phi= ',num2str(Cur(roundn(.5*q,0)))],['\phi= ',num2str(Cur(roundn(.75*q,0)))],['\phi= ',num2str(Cur(roundn(q+1,0)))],'Location','NorthEastOutside');
figure(12)
plot(CUR,EI,'black--','LineWidth',3)
title('# Flextural Rigidity(EI)-Curvature diagram #','Color','b'); 
xlabel('CURVATURE (1/m)');ylabel('EI (kN.m^2)');grid on;
figure(13)
plot(XX,EI,'black','LineWidth',3)
title('# Neuteral axis-Flextural Rigidity(EI) diagram #','Color','b'); 
xlabel('Neuteral axis (mm)');ylabel('EI (kN.m^2)');grid on;
figure(14)
plot(Cur,Mom,X,Y,'r--','LineWidth',3)
title(['# COMPOSITE BEAM SECTION MOMENT-CURVATURE DIAGRAM #','  EI_e_l_a : ',int2str(Elastic_EI),' (kN.m^2)  -  EI_p_l_a : ',int2str(Plastic_EI),' (kN.m^2)'],'Color','b'); 
xlabel('CURVATURE (1/m)');ylabel('MOMENT (kN-m)')
legend('Composite','Composite - bilinear fitted','Location','NorthEastOutside');grid on;
figure(15)
plot(Cur,Mom,CurP,MomP,'r--','LineWidth',3)
title(['# COMPOSITE BEAM CONCRETE SECTION MOMENT-CURVATURE DIAGRAM #'],'Color','b'); 
xlabel('CURVATURE (1/m)');ylabel('MOMENT (kN-m)')
legend('Analysis','SAP2000','Location','NorthEastOutside');grid on;
figure(16)
plot(Cur,Mom,Curvature*1000,Moment*10^-6,'r--','LineWidth',3)
title(['# UNCONFINED CONCRETE SECTION MOMENT-CURVATURE DIAGRAM #'],'Color','b'); 
xlabel('CURVATURE (1/m)');ylabel('MOMENT (kN-m)')
legend('MATLAB','C++','Location','NorthEastOutside');grid on;
%% Output Data to .txt file
fid = fopen('FiberCompositeUnconfinedBeamSectionMomentCurvature-OutPut.txt','w');
fprintf(fid,'*************************************************************\n');
fprintf(fid,'*                  >> IN THE NAME OF GOD <<                 *\n');
fprintf(fid,'* Moment-Curvature Analysis of Composite Unconfined Beam    *\n');
fprintf(fid,'* Section                                                   *\n');
fprintf(fid,'*-----------------------------------------------------------*\n');
fprintf(fid,'*     This program is written by salar delavar ghashghaei   *\n');  
fprintf(fid,'*          E-mail:salar.d.ghashghaei@gmail.com              *\n');
fprintf(fid,'*-----------------------------------------------------------*\n');
fprintf(fid,'*Unit: Newton-Milimeter                                     *\n');
fprintf(fid,'*Given:Section Properties , Steel Section properties ,      *\n');
fprintf(fid,'*Calculate: Moment-Curavture                                *\n');
fprintf(fid,'*************************************************************\n');
fprintf(fid,'*   _    ______________________________________             *\n');
fprintf(fid,'*   |   |                                      |            *\n');
fprintf(fid,'*       |     #     #     #     #    #    #    |            *\n');
fprintf(fid,'*       |     #                           #    |            *\n');
fprintf(fid,'* Beff  |    As1   As2   As3   As4  As5  As6   |            *\n');
fprintf(fid,'*       |     #                           #    |            *\n');
fprintf(fid,'*   |   |     #     #     #     #    #    #    |            *\n');
fprintf(fid,'*   _   |______________________________________|            *\n');
fprintf(fid,'*       |<-               Teff               ->|            *\n');
fprintf(fid,'*       |<-d1->|                                            *\n');
fprintf(fid,'*       |<-  d2   ->|                                       *\n');
fprintf(fid,'*       |<-     d3      ->|                                 *\n');
fprintf(fid,'*       |<-        d4          ->|                          *\n');
fprintf(fid,'*       |<-            d5          ->|                      *\n');
fprintf(fid,'*       |<-               d6             >|                 *\n');
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
fprintf(fid,'*    X                                                      *\n');
fprintf(fid,'*    ^                                                      *\n');
fprintf(fid,'*    |             (Moment - Curvature along X axis)        *\n');
fprintf(fid,'*    |                                                      *\n');
fprintf(fid,'*    +----> Y                                               *\n');
fprintf(fid,'*************************************************************\n');
fprintf(fid,'\n');
fprintf(fid,'+==============================+\n');
fprintf(fid,'=Steel Section Moment-Curvature=\n');
fprintf(fid,'     Curvature    Moment        \n');
fprintf(fid,'       (1/m)      (kN.m)        \n');
fprintf(fid,'--------------------------------\n');
fprintf(fid,'%10.5f %10.3f\n',[Cur;Mom]);
fprintf(fid,'+=============================+\n');
fprintf(fid,'= Steel Section curve fitted  =\n');
fprintf(fid,'     Curvature    Moment       \n');
fprintf(fid,'       (1/m)      (kN.m)       \n');
fprintf(fid,'-------------------------------\n');
fprintf(fid,'%10.4f %10.3f\n',[X;Y]);
fprintf(fid,'+=============================+\n');
fprintf(fid,'\n');
fprintf(fid,'+==============================================================================+\n');
fprintf(fid,' Increment   Top strain   Neuteral axis(x)    Curvature   Flextural Rigidity(EI)\n');
fprintf(fid,'================================================================================\n');
fprintf(fid,'  (i)           (1)           (mm)              (1/m)            (kN.m^2)       \n');
fprintf(fid,'--------------------------------------------------------------------------------\n');
fprintf(fid,'%4.0f %17.5f %13.2f %17.6f %17.2f\n',[I;ES;XX;CUR;EI]);
fprintf(fid,'+===============================================================================+\n');
fclose(fid);