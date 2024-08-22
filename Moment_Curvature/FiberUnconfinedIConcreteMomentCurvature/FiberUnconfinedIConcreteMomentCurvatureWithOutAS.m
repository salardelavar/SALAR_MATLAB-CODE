%***********************************************************%
%                  >> IN THE NAME OF GOD <<                 %
% Moment-Curvature Analysis of Unconfined I concrete section%
% With and Without Steel Reinforcement                      %
%-----------------------------------------------------------%
%     This program is written by salar delavar ghashghaei   %  
%          E-mail:salar.d.ghashghaei@gmail.com              %
%-----------------------------------------------------------%
%Unit: Newton-Milimeter                                     %
%Given:Section Properties , Concrete properties ,           %
% Reinforcing steel properties                              %
%Calculate: Moment-Curavture                                %
% Note: No limit for accounting plurality steel rebar       %
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
%        ______________________________________             %
%       |                                      |            %
%       |     #     #     #     #    #    #    |            %
%       |     #                           #    |            %
%       |    As1   As2   As3   As4  As5  As6   |            %
%       |     #                           #    |            %
%       |     #     #     #     #    #    #    |            %
%       |______________________________________|            %
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
clear all;close all;clc;
%% Section Properties
tf1=150;% [mm] I section thickness on Top flange
bf1=1100;% [mm] I section width on Top flange
tw=200;% [mm] I section thickness of Web
hw=1900;% [mm] Height of web
tf2=150;% [mm] I section thickness on Bottom flange
bf2=1100;% [mm] I section width on Bottom flange
h=tf1+tf2+hw;% [mm] Height of Section
%As:  As1      As2     As3     As4    As5      As6
As=[3500 1000 1000 1000 1000 3500]; 
%d:d1  d2  d3  d4  d5  d6 
d=[75 625 1050 1625 2050 2125];
%% Concrete Properties
fc =25;% [N/mm^2] Unconfined concrete strength 
ecu=0.004;% Ultimate concrete strain
Ec=5000*sqrt(fc);
ec0=(2*fc)/Ec;
fct=-0.7*sqrt(fc);% Concrete tension stress
ect1=(2*fct)/Ec;ect2=(2.625*fct)/Ec;ect3=(9.292*fct)/Ec;% Concrete tension strain
%% Reinforcing steel Properties
fy =400;% [N/mm^2] Yield strength of reinforcing steel
Es =2e5;% [N/mm^2] Modulus of elasticity of steel
fu=1.5*fy;% Ultimate steel stress
ey=fy/Es;% Yield steel strain
esh=0.01;% Strain at steel strain-hardening
esu=0.09;% Ultimate steel strain
Esh=(fu-fy)/(esu-esh);
N=h*2;% Number of concrete Fiber
itermax = 4000;% maximum number of iterations
tolerance = 1e-12; % specified tolerance for convergence
x=.5*h;% initial guess of Neuteral axis
%%% monitor cpu time
starttime = cputime;
%% ------------------ Newton Method Procedure - With Reinforcement ------------------------%
disp('################################################');
disp('# Moment-Curvature Analysis With Reinforcement #');
disp('################################################');
An=size(As,2);
N1=N/10;% Number of steel section top flange Fiber  
N2=8*N/10;% Number of steel section  web Fiber
N3=N/10;% Number of steel section bottom flange Fiber 
R1=(1/N1);R2=(1/N2);R3=(1/N3);
for k=1:N1;c1(k)=(.5*R1+R1*(k-1))*tf1;end % distance of each top flange Fiber 
for k=1:N2;c2(k)=tf1+(.5*R2+R2*(k-1))*hw;end % distance of each web Fiber
for k=1:N3;c3(k)=tf1+hw+(.5*R3+R3*(k-1))*tf2;end % distance of each bottom flange Fiber
c=[c1 c2 c3];
EC=[abs(ect1) abs(ect2) .33*abs(ect3) .67*abs(ect3) .8*abs(ect3) .9*abs(ect3) abs(ect3) .4*ecu .5*ecu .6*ecu .7*ecu .8*ecu .9*ecu ecu];q=size(EC,2);
% EC=abs(ect1):0.00005:ecu;q=size(EC,2);
for j=1:q;
    eC=EC(j);
    it = 0; % initialize iteration count
    residual = 100; % initialize residual
    while (residual > tolerance)
       for u=1:An;
       es=eC*(x-d(u))/x;
    %---------------- As -------------%
       if and(es>0,es<ey)
        fs(u)=Es*es;
        fstan(u)=(Es*eC*d(u))/(x)^2;
       elseif and(es<0,es>(-ey))
        fs(u)=Es*es;
        fstan(u)=(Es*eC*d(u))/(x)^2;
       elseif  and(es>=ey,es<esh)
        fs(u)=fy;
        fstan(u)=0;
       elseif and(es<=(-ey),es>(-esh))
        fs(u)=-fy;
        fstan(u)=0;
        elseif  and(es>=esh,es<esu)
        fs(u)=fy+Esh*(abs(es)-esh);
        fstan(u)=(Esh*eC*d(u))/(x)^2;
       elseif and(es<=(-esh),es>(-esu))
        fs(u)=-fy-Esh*(abs(es)-esh);
        fstan(u)=(Esh*eC*d(u))/(x)^2;
       elseif or(es>=esu,es<=(-esu))
        fs(u)=0;
        fstan(u)=0;
     end
        Fs(u)=As(u)*fs(u);
        Fstan(u)=As(u)*fstan(u);% tangent steel force
        if and(es>0,es<ec0)% in this step: Unconfined concrete force in rebar area is omitted (F=As*fc)
        Cs=fc*((2*es/ec0)-(es/ec0)^2);
        Ctans=((2*fc)/(ec0^2*x^3))*(ec0*eC*d(u)*x-d(u)*eC^2*x+2*eC^2*d(u)^2);
      elseif and(es>=ec0,es<ecu)
        Cs=fc*(1-(0.15*(es-ec0)/(ecu-ec0)));
        Ctans=-(3*eC*d(u)*fc)/(20*(ecu-ec0)*x^2);
      elseif es >= ecu
         Cs=0; 
         Ctans=0;
      elseif and(es<0,es>=ect1)
        Cs=-0.5*Ec*es;
        Ctans=-(0.5*Ec*d(u)*eC)/x^2;
      elseif and(es<ect1,es>=ect2)
        Cs=fct+(0.5*fct/(ect2-ect1))*(es-ect1);
        Ctans=+(0.5*fct*eC*d(u))/((ect2-ect1)*x^2);
      elseif and(es<ect2,es>=ect3)
        Cs=-(.5*fct-(0.5*fct/(ect3-ect2))*(es-ect2));
        Ctans=+(0.5*fct*eC*d(u))/((ect3-ect2)*x^2);
      elseif es<ect3
        Cs=0;
        Ctans=0;
        end
        CS(u)=Cs;CtanS(u)=Ctans;
        Fss(u)=-As(u)*CS(u);
        Fstans(u)=-As(u)*CtanS(u);% tangent Minus of concrete force
    end                
    for z=1:N;% in this step: concrete force for each fiber is calculated
     %-------------- Cc --------------%
     ec=eC*(x-c(z))/x;
     if and(z>=1,z<=roundn(N/10,0)) % TOP FLANGE
     if and(ec>0,ec<ec0)
        C=fc*((2*ec/ec0)-(ec/ec0)^2);
        Ctan=((2*fc)/(ec0^2*x^3))*(ec0*eC*c(z)*x-c(z)*eC^2*x+2*eC^2*c(z)^2);
      elseif and(ec>=ec0,ec<ecu)
        C=fc*(1-(0.15*(ec-ec0)/(ecu-ec0)));
        Ctan=-(3*eC*c(z)*fc)/(20*(ecu-ec0)*x^2);
      elseif ec >= ecu
         C=0; 
         Ctan=0;
      elseif and(ec<0,ec>=ect1)
        C=0.5*Ec*ec;
        Ctan=(0.5*Ec*c(z)*eC)/x^2;
      elseif and(ec<ect1,ec>=ect2)
        C=fct-(0.5*fct/(ect2-ect1))*(ec-ect1);
        Ctan=-(0.5*fct*eC*c(z))/((ect2-ect1)*x^2);
      elseif and(ec<ect2,ec>=ect3)
        C=.5*fct-(0.5*fct/(ect3-ect2))*(ec-ect2);
        Ctan=-(0.5*fct*eC*c(z))/((ect3-ect2)*x^2);
      elseif ec<ect3
        C=0;
        Ctan=0;
     end
     Cc(z)=bf1*R1*tf1*C;
     Cctan(z)=bf1*R1*tf1*Ctan;% tangent concrete force
     elseif and(z>roundn(N/10,0),z<=roundn(9*N/10,0)) % MIDDLE WEB
     if and(ec>0,ec<ec0)
        C=fc*((2*ec/ec0)-(ec/ec0)^2);
        Ctan=((2*fc)/(ec0^2*x^3))*(ec0*eC*c(z)*x-c(z)*eC^2*x+2*eC^2*c(z)^2);
      elseif and(ec>=ec0,ec<ecu)
        C=fc*(1-(0.15*(ec-ec0)/(ecu-ec0)));
        Ctan=-(3*eC*c(z)*fc)/(20*(ecu-ec0)*x^2);
      elseif ec >= ecu
         C=0; 
         Ctan=0;
      elseif and(ec<0,ec>=ect1)
        C=0.5*Ec*ec;
        Ctan=(0.5*Ec*c(z)*eC)/x^2;
      elseif and(ec<ect1,ec>=ect2)
        C=fct-(0.5*fct/(ect2-ect1))*(ec-ect1);
        Ctan=-(0.5*fct*eC*c(z))/((ect2-ect1)*x^2);
      elseif and(ec<ect2,ec>=ect3)
        C=.5*fct-(0.5*fct/(ect3-ect2))*(ec-ect2);
        Ctan=-(0.5*fct*eC*c(z))/((ect3-ect2)*x^2);
      elseif ec<ect3
        C=0;
        Ctan=0;
     end
     Cc(z)=tw*R2*hw*C;
     Cctan(z)=tw*R2*hw*Ctan;% tangent concrete force
     elseif and(z>roundn(9*N/10,0),z<=N) % BELOW FLANGE
              if and(ec>0,ec<ec0)
        C=fc*((2*ec/ec0)-(ec/ec0)^2);
        Ctan=((2*fc)/(ec0^2*x^3))*(ec0*eC*c(z)*x-c(z)*eC^2*x+2*eC^2*c(z)^2);
      elseif and(ec>=ec0,ec<ecu)
        C=fc*(1-(0.15*(ec-ec0)/(ecu-ec0)));
        Ctan=-(3*eC*c(z)*fc)/(20*(ecu-ec0)*x^2);
      elseif ec >= ecu
         C=0; 
         Ctan=0;
      elseif and(ec<0,ec>=ect1)
        C=0.5*Ec*ec;
        Ctan=(0.5*Ec*c(z)*eC)/x^2;
      elseif and(ec<ect1,ec>=ect2)
        C=fct-(0.5*fct/(ect2-ect1))*(ec-ect1);
        Ctan=-(0.5*fct*eC*c(z))/((ect2-ect1)*x^2);
      elseif and(ec<ect2,ec>=ect3)
        C=.5*fct-(0.5*fct/(ect3-ect2))*(ec-ect2);
        Ctan=-(0.5*fct*eC*c(z))/((ect3-ect2)*x^2);
      elseif ec<ect3
        C=0;
        Ctan=0;
     end
     Cc(z)=bf2*R3*tf2*C;
     Cctan(z)=bf2*R3*tf2*Ctan;% tangent concrete force
     end
     Ss(z)=ec;SS(j,z)=Ss(z);% Concrete Fiber Strain
     CFS(j,z)=C;% Concrete Fiber Stress
    end                
    %----------------------------------%
    FsTOTAL=sum(Fs);CcTOTAL=sum(Cc);FssTOTAL=sum(Fss);A=CcTOTAL+FsTOTAL;
    FsTOTAL_tan=sum(Fstan);CcTOTAL_tan=sum(Cctan);FssTOTAL_tan=sum(Fstans);A_tan=CcTOTAL_tan+FsTOTAL_tan;
    dx = A_tan^-1 *(-A);
    residual = max(abs(dx)); % evaluate residual
        it = it + 1; % increment iteration count
        x = x+dx; % update x
        if it == itermax % stop the the analysis of this step please of Convergence
          fprintf('(-)Increment %1.0f : trail iteration reached to Ultimate %1.0f - strain: %1.6f - error: [%1.2f]\n',j,it,eC,A)
             disp('    ## The solution for this step is not converged. Please check your model ##') 
            break
        end
    end
    if it == itermax;break;end % stop the analysis at all because last Convergence
    for u=1:An;e(u)=x-d(u);end % distance of each rebar from Neuteral axis
     for k=1:N;cc(k)=x-c(k);end %distance of each concrete fiber from Neuteral axis
        if it < itermax% iteration control
        fprintf('(+)Increment %1.0f : It is converged in %1.0f iterations - strain: %1.6f - x: %1.2f - Phi: %1.5f - Moment: %1.2f\n',j,it,eC,x,(eC/x)*1000,(Fs*e'+Cc*cc')*10^-6)
        end
        TSCSs(j)=eC*(x-d(1))/x;TSCS(j)=fs(1);% Top Steel compression strain-stress
        BSCSs(j)=eC*(x-d(6))/x;BSCS(j)=fs(6);% Bottom Steel compression strain-stress
if and(eC>0,eC<ec0);ECCC(j)=eC;TUCS(j)=fc*((2*eC/ec0)-(eC/ec0)^2);elseif and(eC>=ec0,eC<ecu);ECCC(j)=eC;TUCS(j)=fc*(1-(0.15*(eC-ec0)/(ecu-ec0)));end% Top Unconfined compression stress 
% Calculate Moment and Curavture
Cur(j)=(eC/x)*1000;XX(j)=x;AA(j)=A;
Mom(j)=(Fs*e'+Fss*e'+Cc*cc')*10^-6;
CUR(j)=Cur(j);I(j)=j;IT(j)=it;DX(j)=dx;
end
Cur=[0 Cur];Mom=[0 Mom];
s=size(Cur,2);for i=1:s-1;EI(i)=(Mom(i+1)-Mom(i))/(Cur(i+1)-Cur(i));end % Flextural Rigidity
if roundn(eC,-5) == ecu;fprintf('\n      ## Unconfined Concrete Strain Reached to Ultimate Strain: %1.4f ## \n\n',eC);end
%% UnConfined bilinear fitting
SIZE=size(Mom,2);
for i=1:SIZE-1;
    hh(i) = Cur(i+1)-Cur(i);
    Aa(i)=(Mom(i)+Mom(i+1))*0.5*hh(i);
end
Area=sum(Aa);k0 =Mom(2)/Cur(2);
fiy = (Mom(i+1)*max(Cur)*0.5-Area)/(Mom(i+1)*0.5 - k0*max(Cur)*0.5);
My = k0*fiy;
X = [0 fiy max(Cur)];Y = [0 My Mom(i+1)];
disp('+=================================+');
disp('= Unconfined curve fitted With As =');
disp('  Curvature          Moment        ');
disp('    (1/m)            (kN.m)        ');
disp('-----------------------------------');
disp([X' Y']);
disp('+=================================+');
%% ------------------ Newton Method Procedure - Without Reinforcement ------------------------%
disp('####################################################');
disp('#  Moment-Curvature Analysis Without Reinforcement #');
disp('####################################################');
R=(1/N);for k=1:N;co(k)=(.5*R+R*(k-1))*h;end
% ECo=[.2*abs(ect1) .4*abs(ect1) .6*abs(ect1) .8*abs(ect1) abs(ect1) .2*abs(ect2) .4*abs(ect2) .6*abs(ect2) .8*abs(ect2) abs(ect2) .2*abs(ect3) .4*abs(ect3) .6*abs(ect3) .8*abs(ect3) abs(ect3)-(1^-5)];qo=size(ECo,2);
ECo=10^-10:5*10^-6:0.000362;qo=size(ECo,2);
for jj=1:qo;
    eCo=ECo(jj);
    it = 0; % initialize iteration count
    residual = 100; % initialize residual
    while (residual > tolerance)                
    for z=1:N;% in this step: concrete force for each fiber is calculated
     %-------------- Cc --------------%
     eco=eCo*(x-co(z))/x;
     if and(z>=1,z<=roundn(N/10,0)) % TOP FLANGE
     if and(eco>0,eco<ec0)
        Co=fc*((2*eco/ec0)-(eco/ec0)^2);
        Ctano=((2*fc)/(ec0^2*x^3))*(ec0*eCo*co(z)*x-co(z)*eCo^2*x+2*eCo^2*co(z)^2);
      elseif and(eco>=ec0,eco<ecu)
        Co=fc*(1-(0.15*(eco-ec0)/(ecu-ec0)));
        Ctano=-(3*eCo*co(z)*fc)/(20*(ecu-ec0)*x^2);
      elseif eco >= ecu
         Co=0; 
         Ctano=0;
      elseif and(eco<0,eco>=ect1)
        Co=0.5*Ec*eco;
        Ctano=(0.5*Ec*co(z)*eCo)/x^2;
      elseif and(eco<ect1,eco>=ect2)
        Co=fct-(0.5*fct/(ect2-ect1))*(eco-ect1);
        Ctano=-(0.5*fct*eCo*co(z))/((ect2-ect1)*x^2);
      elseif and(eco<ect2,eco>=ect3)
        Co=.5*fct-(0.5*fct/(ect3-ect2))*(eco-ect2);
        Ctano=-(0.5*fct*eCo*co(z))/((ect3-ect2)*x^2);
      elseif eco<ect3
        Co=0;
        Ctano=0;
     end
     Cco(z)=bf1*R1*tf1*Co;
     Cctano(z)=bf1*R1*tf1*Ctano;% tangent concrete force
     elseif and(z>roundn(N/10,0),z<=roundn(9*N/10,0)) % MIDDLE WEB
          if and(eco>0,eco<ec0)
        Co=fc*((2*eco/ec0)-(eco/ec0)^2);
        Ctano=((2*fc)/(ec0^2*x^3))*(ec0*eCo*co(z)*x-co(z)*eCo^2*x+2*eCo^2*co(z)^2);
      elseif and(eco>=ec0,eco<ecu)
        Co=fc*(1-(0.15*(eco-ec0)/(ecu-ec0)));
        Ctano=-(3*eCo*co(z)*fc)/(20*(ecu-ec0)*x^2);
      elseif eco >= ecu
         Co=0; 
         Ctano=0;
      elseif and(eco<0,eco>=ect1)
        Co=0.5*Ec*eco;
        Ctano=(0.5*Ec*co(z)*eCo)/x^2;
      elseif and(eco<ect1,eco>=ect2)
        Co=fct-(0.5*fct/(ect2-ect1))*(eco-ect1);
        Ctano=-(0.5*fct*eCo*co(z))/((ect2-ect1)*x^2);
      elseif and(eco<ect2,eco>=ect3)
        Co=.5*fct-(0.5*fct/(ect3-ect2))*(eco-ect2);
        Ctano=-(0.5*fct*eCo*co(z))/((ect3-ect2)*x^2);
      elseif eco<ect3
        Co=0;
        Ctano=0;
     end
     Cco(z)=tw*R2*hw*Co;
     Cctano(z)=tw*R2*hw*Ctano;% tangent concrete force
     elseif and(z>roundn(9*N/10,0),z<=N) % BELOW FLANGE
              if and(eco>0,eco<ec0)
        Co=fc*((2*eco/ec0)-(eco/ec0)^2);
        Ctano=((2*fc)/(ec0^2*x^3))*(ec0*eCo*co(z)*x-co(z)*eCo^2*x+2*eCo^2*co(z)^2);
      elseif and(eco>=ec0,eco<ecu)
        Co=fc*(1-(0.15*(eco-ec0)/(ecu-ec0)));
        Ctano=-(3*eCo*co(z)*fc)/(20*(ecu-ec0)*x^2);
      elseif eco >= ecu
         Co=0; 
         Ctano=0;
      elseif and(eco<0,eco>=ect1)
        Co=0.5*Ec*eco;
        Ctano=(0.5*Ec*co(z)*eCo)/x^2;
      elseif and(eco<ect1,eco>=ect2)
        Co=fct-(0.5*fct/(ect2-ect1))*(eco-ect1);
        Ctano=-(0.5*fct*eCo*co(z))/((ect2-ect1)*x^2);
      elseif and(eco<ect2,eco>=ect3)
        Co=.5*fct-(0.5*fct/(ect3-ect2))*(eco-ect2);
        Ctano=-(0.5*fct*eCo*co(z))/((ect3-ect2)*x^2);
      elseif eco<ect3
        Co=0;
        Ctano=0;
     end
     Cco(z)=bf2*R3*tf2*Co;
     Cctano(z)=bf2*R3*tf2*Ctano;% tangent concrete force
     end
     Sso(z)=eco;SSo(jj,z)=Sso(z);% Concrete Fiber Strain
     CFSo(jj,z)=Co;% Concrete Fiber Stress
    end                
    %----------------------------------%
    CcTOTAL=sum(Cco);A=CcTOTAL;
    CcTOTAL_tan=sum(Cctano);A_tan=CcTOTAL_tan;
    dx = A_tan^-1 *(-A);
    residual = max(abs(dx)); % evaluate residual
        it = it + 1; % increment iteration count
        x = x+dx; % update x
        if it == itermax % stop the the analysis of this step please of Convergence
          fprintf('(-)Increment %1.0f : trail iteration reached to Ultimate %1.0f - strain: %1.6f - error: [%1.2f]\n',jj,it,eCo,A)
             disp('    ## The solution for this step is not converged. Please check your model ##') 
            %break
        end
    end
    if it == itermax;break;end % stop the analysis at all because last Convergence     
     for k=1:N;cco(k)=x-co(k);end %distance of each concrete fiber from Neuteral axis
        if it < itermax% iteration control
        fprintf('(+)Increment %1.0f : It is converged in %1.0f iterations - strain: %1.6f - x: %1.2f - Phi: %1.5f - Moment: %1.2f\n',jj,it,eCo,x,(eCo/x)*1000,(Cco*cco')*10^-6)
        end
if and(eCo>0,eCo<ec0);ECCCo(jj)=eCo;TUCSo(jj)=fc*((2*eCo/ec0)-(eCo/ec0)^2);elseif and(eCo>=ec0,eCo<ecu);ECCCo(jj)=eCo;TUCSo(jj)=fc*(1-(0.15*(eCo-ec0)/(ecu-ec0)));end% Top Unconfined compression stress 
% Calculate Moment and Curavture
Curo(jj)=(eCo/x)*1000;XXo(jj)=x;AAo(jj)=A;
Momo(jj)=(Cco*cco')*10^-6;
CURo(jj)=Curo(jj);Io(jj)=jj;ITo(jj)=it;DXo(jj)=dx;
end
Curo=[0 Curo];Momo=[0 Momo];
so=size(Curo,2);for i=1:so-1;EIo(i)=(Momo(i+1)-Momo(i))/(Curo(i+1)-Curo(i));end % Flextural Rigidity
if roundn(eCo,-5) == abs(ect3);fprintf('\n      ## Without Reinforcement Unconfined Concrete Strain Reached to Ultimate Tension Strain: %1.4f ## \n\n',eCo);end
%% UnConfined bilinear fitting - Without Reinforcement
SIZEo=size(Momo,2);
for i=1:SIZEo-1;
    hho(i) = Curo(i+1)-Curo(i);
    Aao(i)=(Momo(i)+Momo(i+1))*0.5*hho(i);
end
Areao=sum(Aao);k0o =Momo(2)/Curo(2);
fiyo = (Momo(i+1)*max(Curo)*0.5-Areao)/(Momo(i+1)*0.5 - k0o*max(Curo)*0.5);
Myo = k0o*fiyo;
Xo = [0 fiyo max(Curo)];Yo = [0 Myo Momo(i+1)];
disp('+====================================+');
disp('= Unconfined curve fitted Without As =');
disp('  Curvature             Moment        ');
disp('    (1/m)               (kN.m)        ');
disp('--------------------------------------');
disp([Xo' Yo']);
disp('+=====================================+');
%% EI and Ductility_Rito of  Unconfined Section 
Elastic_EI=Y(2)/X(2);
Plastic_EI=(Y(3)-Y(2))/(X(3)-X(2));
Ductility_Rito=X(3)/X(2);
Elastic_EIo=Yo(2)/Xo(2);
Plastic_EIo=(Yo(3)-Yo(2))/(Xo(3)-Xo(2));
Ductility_Ritoo=Xo(3)/Xo(2);
fprintf('+-----------------------------------------------------------------+\n')
fprintf(' With Reinforcement - Elastic EI : %5.2f (kN.m^2)\n',Elastic_EI)
fprintf(' With Reinforcement -Plastic EI : %5.2f (kN.m^2)\n',Plastic_EI)
fprintf(' With Reinforcement - Unconfined Section Ductility Rito : %5.2f\n',Ductility_Rito)
fprintf(' Without Reinforcement - Elastic EI : %5.2f (kN.m^2)\n',Elastic_EIo)
fprintf(' Without Reinforcement -Plastic EI : %5.2f (kN.m^2)\n',Plastic_EIo)
fprintf(' Without Reinforcement - Unconfined Section Ductility Rito : %5.2f\n',Ductility_Ritoo)
fprintf('+-----------------------------------------------------------------+\n')
%%%  print time of computation
totaltime = cputime - starttime;
fprintf('\nTotal time (s): %7.4f \n\n',totaltime)
%% imaging
figure(1)
IMAGE1=imread('FiberUnconfinedConcreteMomentCurvature-image1.jpg');
image(IMAGE1);axis image;axis off;
figure(2)
IMAGE2=imread('FiberUnconfinedIConcreteMomentCurvature.jpg');
image(IMAGE2);axis image;axis off;
%% Plot
figure(3)
p1=plot(I,DX,Io,DXo,'r--');grid on;set(p1,'LineWidth',2);
xlabel('increment');ylabel('Residual');
legend('With Reinf.','Without Reinf.','Location','NorthEastOutside');grid on;
title('Residual-increment diagram','color','b');
figure(4)
p1=plot(I,IT,Io,ITo,'r--');grid on;set(p1,'LineWidth',2);
xlabel('increment');ylabel('Iteration');
legend('With Reinf.','Without Reinf.','Location','NorthEastOutside');grid on;
title('Iteration-increment diagram','color','b');
figure(5)
plot(EC,XX,ECo,XXo,'r--','LineWidth',2);xlabel('Top concrete strain ');ylabel('Neuteral axis (mm)')
legend('With Reinf.','Without Reinf.','Location','NorthEastOutside');grid on;
title(' Plotting of the Top concrete strain and Neuteral axis ','Color','b');grid on
figure(6)
plot(EC,CUR,ECo,CURo,'r--','LineWidth',2);xlabel('Top concrete strain ');ylabel('Curvature (1/m)')
legend('With Reinf.','Without Reinf.','Location','NorthEastOutside');grid on;
title(' Plotting of the Top concrete strain and Curvature ','Color','b');grid on
figure(7)
plot(ECCC,TUCS,ECCCo,TUCSo,'r--','LineWidth',3);xlabel('Top fiber concrete strain ');ylabel('Top fiber concrete stress (N/mm^2)')
legend('With Reinf.','Without Reinf.','Location','NorthEastOutside');grid on;
title(' Plotting of the Top fiber unconfined concrete strain-stress during the analysis','Color','b');grid on
figure(8)
plot(SS(1,:),(h-c),SS(roundn(.25*q,0),:),(h-c),SS(roundn(.5*q,0),:),(h-c),SS(roundn(.75*q,0),:),(h-c),SS(roundn(q,0),:),(h-c),'g--','LineWidth',3);xlabel('Concrete strain (mm/mm)');ylabel('Height (mm)')
title(' Plotting of the concrete fiber strain and Heigh of section - Negative value: Tension - Positive value: Compression','Color','b');grid on
legend(['\phi= ',num2str(Cur(2))],['\phi= ',num2str(Cur(roundn(.25*q,0)))],['\phi= ',num2str(Cur(roundn(.5*q,0)))],['\phi= ',num2str(Cur(roundn(.75*q,0)))],['\phi= ',num2str(Cur(roundn(q+1,0)))],'Location','NorthEastOutside');
figure(9)
plot(SSo(1,:),(h-co),SSo(roundn(.25*qo,0),:),(h-co),SSo(roundn(.5*qo,0),:),(h-co),SSo(roundn(.75*qo,0),:),(h-co),SSo(roundn(qo,0),:),(h-co),'g--','LineWidth',3);xlabel('Concrete strain (mm/mm)');ylabel('Height (mm)')
title(' Plotting of the concrete fiber strain and Heigh of section - Negative value: Tension - Positive value: Compression','Color','b');grid on
legend(['\phi= ',num2str(Curo(2))],['\phi= ',num2str(Curo(roundn(.25*qo,0)))],['\phi= ',num2str(Curo(roundn(.5*qo,0)))],['\phi= ',num2str(Curo(roundn(.75*qo,0)))],['\phi= ',num2str(Curo(roundn(qo+1,0)))],'Location','NorthEastOutside');
figure(10)
plot(CFS(1,:),(h-c),CFS(roundn(.25*q,0),:),(h-c),CFS(roundn(.5*q,0),:),(h-c),CFS(roundn(.75*q,0),:),(h-c),CFS(roundn(q,0),:),(h-c),'g--','LineWidth',3);xlabel('Concrete stress (Mpa)');ylabel('Height (mm)')
title(' Plotting of the concrete fiber stress and Heigh of section - Negative value: Tension - Positive value: Compression','Color','b');grid on
legend(['\phi= ',num2str(Cur(2))],['\phi= ',num2str(Cur(roundn(.25*q,0)))],['\phi= ',num2str(Cur(roundn(.5*q,0)))],['\phi= ',num2str(Cur(roundn(.75*q,0)))],['\phi= ',num2str(Cur(roundn(q+1,0)))],'Location','NorthEastOutside');
figure(11)
plot(CFSo(1,:),(h-co),CFSo(roundn(.25*qo,0),:),(h-co),CFSo(roundn(.5*qo,0),:),(h-co),CFSo(roundn(.75*qo,0),:),(h-co),CFSo(roundn(qo,0),:),(h-co),'g--','LineWidth',3);xlabel('Concrete stress (Mpa)');ylabel('Height (mm)')
title(' Plotting of the concrete fiber stress and Heigh of section - Negative value: Tension - Positive value: Compression','Color','b');grid on
legend(['\phi= ',num2str(Curo(2))],['\phi= ',num2str(Curo(roundn(.25*qo,0)))],['\phi= ',num2str(Curo(roundn(.5*qo,0)))],['\phi= ',num2str(Curo(roundn(.75*qo,0)))],['\phi= ',num2str(Curo(roundn(qo+1,0)))],'Location','NorthEastOutside');
figure(12)
plot(TSCSs,TSCS,BSCSs,BSCS,'r--','LineWidth',3);xlabel('Steel strain (\epsilon_s)');ylabel('Steel stress (N/mm^2)')
title(' Plotting of the Top  and Bottom steel strain-stress during the analysis','Color','b');
legend('Top steel','Bottom steel','Location','NorthEastOutside');grid on;
figure(13)
plot(CUR,EI,CURo,EIo,'r--','LineWidth',3)
title('# UNCONFINED Flextural Rigidity(EI)-CURVATURE DIAGRAM #','Color','b');
legend('With Reinf.','Without Reinf.','Location','NorthEastOutside');grid on;
xlabel('CURVATURE (1/m)');ylabel('EI (kN.m^2)');grid on;
figure(12)
plot(XX,EI,XXo,EIo,'r--','LineWidth',3)
title('# Neuteral Axis-Flextural Rigidity(EI) diagram #','Color','b');
legend('With Reinf.','Without Reinf.','Location','NorthEastOutside');grid on;
xlabel('Neuteral axis (mm)');ylabel('EI (kN.m^2)');grid on;
figure(13)
plot(Cur,Mom,X,Y,'r--','LineWidth',3)
title(['# WITH REINFORCEMENT UNCONFINED MOMENT CURVATURE DIAGRAM #',' -  EI_e_l_a : ',int2str(Elastic_EI),' (kN.m^2)  -  EI_p_l_a : ',int2str(Plastic_EI),' (kN.m^2)'],'Color','b'); 
xlabel('CURVATURE (1/m)');ylabel('MOMENT (kN-m)')
legend('Unconfined','Unconfined bilinear fitted','Location','NorthEastOutside');grid on;
figure(14)
plot(Curo,Momo,Xo,Yo,'r--','LineWidth',3)
title(['# WITHIOUT REINFORCEMENT UNCONFINED MOMENT CURVATURE DIAGRAM #',' -  EI_e_l_a : ',int2str(Elastic_EIo),' (kN.m^2)  -  EI_p_l_a : ',int2str(Plastic_EIo),' (kN.m^2)'],'Color','b'); 
xlabel('CURVATURE (1/m)');ylabel('MOMENT (kN-m)')
legend('Unconfined','Unconfined bilinear fitted','Location','NorthEastOutside');grid on;
figure(15)
plot(Cur,Mom,Curo,Momo,'r--','LineWidth',3)
title('# WITH AND WITHIOUT REINFORCEMENT UNCONFINED MOMENT CURVATURE DIAGRAM #','Color','b'); 
xlabel('CURVATURE (1/m)');ylabel('MOMENT (kN-m)')
legend('With Reinf.','Without Reinf.','Location','NorthEastOutside');grid on;
%% Output Data to .txt file
fid = fopen('FiberUnconfinedIConcreteMomentCurvatureWithOutAS-OutPut.txt','w');
fprintf(fid,'*************************************************************\n');
fprintf(fid,'*                  >> IN THE NAME OF GOD <<                 *\n');
fprintf(fid,'*Analysis of Unconfined I concrete section Moment-Curvature *\n');
fprintf(fid,'* With and Without Steel Reinforcement                      *\n');
fprintf(fid,'*-----------------------------------------------------------*\n');
fprintf(fid,'*     This program is written by salar delavar ghashghaei   *\n');  
fprintf(fid,'*          E-mail:salar.d.ghashghaei@gmail.com              *\n');
fprintf(fid,'*-----------------------------------------------------------*\n');
fprintf(fid,'*Unit: Newton-Milimeter                                     *\n');
fprintf(fid,'*Given:Section Properties , Concrete properties ,           *\n');
fprintf(fid,'* Reinforcing steel properties                              *\n');
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
fprintf(fid,'*        ______________________________________             *\n');
fprintf(fid,'*       |                                      |            *\n');
fprintf(fid,'*       |     #     #     #     #    #    #    |            *\n');
fprintf(fid,'*       |     #                           #    |            *\n');
fprintf(fid,'*       |    As1   As2   As3   As4  As5  As6   |            *\n');
fprintf(fid,'*       |     #                           #    |            *\n');
fprintf(fid,'*       |     #     #     #     #    #    #    |            *\n');
fprintf(fid,'*       |______________________________________|            *\n');
fprintf(fid,'*       |<-d1->|                                            *\n');
fprintf(fid,'*       |<-  d2   ->|                                       *\n');
fprintf(fid,'*       |<-     d3      ->|                                 *\n');
fprintf(fid,'*       |<-        d4          ->|                          *\n');
fprintf(fid,'*       |<-            d5          ->|                      *\n');
fprintf(fid,'*       |<-               d6             >|                 *\n');
fprintf(fid,'*************************************************************\n');
fprintf(fid,'\n');
fprintf(fid,'+==========================+\n');
fprintf(fid,'=   Unconfined with As     =\n');
fprintf(fid,'    Curvature    Moment     \n');
fprintf(fid,'      (1/m)      (kN.m)     \n');
fprintf(fid,'----------------------------\n');
fprintf(fid,'%10.5f %10.3f\n',[Cur;Mom]);
fprintf(fid,'+===========================+\n');
fprintf(fid,'\n');
fprintf(fid,'+===========================+\n');
fprintf(fid,'=   Unconfined without As   =\n');
fprintf(fid,'  Curvature    Moment        \n');
fprintf(fid,'    (1/m)      (kN.m)        \n');
fprintf(fid,'-----------------------------\n');
fprintf(fid,'%10.5f %10.3f\n',[Curo;Momo]);
fprintf(fid,'+===========================+\n');
fprintf(fid,'+==========================================================================+\n');
fprintf(fid,'= Unconfined curve fitted  with   As == Unconfined curve fitted without As =\n');
fprintf(fid,'  Curvature             Moment          Curvature             Moment        \n');
fprintf(fid,'    (1/m)               (kN.m)            (1/m)               (kN.m)        \n');
fprintf(fid,'----------------------------------------------------------------------------\n');
fprintf(fid,'%10.5f %10.3f %10.5f %10.3f\n',[X;Y;Xo;Yo]);
fprintf(fid,'+==========================================================================+\n');
fprintf(fid,'\n');
fprintf(fid,'+============================== With Reinforcement ============================+\n');
fprintf(fid,' Increment   Top strain   Neuteral axis(x)    Curvature   Flextural Rigidity(EI)\n');
fprintf(fid,'================================================================================\n');
fprintf(fid,'  (i)           (1)           (mm)              (1/m)            (kN.m^2)       \n');
fprintf(fid,'--------------------------------------------------------------------------------\n');
fprintf(fid,'%4.0f %17.5f %13.2f %17.6f %17.2f\n',[I;EC;XX;CUR;EI]);
fprintf(fid,'+===============================================================================+\n');
fprintf(fid,'+============================ Without Reinforcement ============================+\n');
fprintf(fid,' Increment   Top strain   Neuteral axis(x)    Curvature   Flextural Rigidity(EI)\n');
fprintf(fid,'================================================================================\n');
fprintf(fid,'  (i)           (1)           (mm)              (1/m)            (kN.m^2)       \n');
fprintf(fid,'--------------------------------------------------------------------------------\n');
fprintf(fid,'%4.0f %17.5f %13.2f %17.6f %17.2f\n',[Io;ECo;XXo;CURo;EIo]);
fprintf(fid,'+===============================================================================+\n');
fclose(fid);