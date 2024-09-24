#include <iostream> //std::cout
#include <conio.h>
#include <fstream>
#include <cstdlib>
#include <cmath> //std::pow - sqrt
#include <ConsoleColor.h>
#include <iomanip> // std:: round number
#include <ctime>
using namespace std;
main()
  {
//system("color 8E");	
//system("color 1F");
//system("color 0A");
	double Force1,Force2,Dini,tolerance,residual;
	double DATA[184000],CHECK[50];
	int m,i,j,it,itermax,N,M,ZA,X=0;	
    ifstream INo;INo.open("PushoverAnalysisNonlinear8SpringsDC2dof-input.txt");
    INo>>Force1>>Force2>>Dini>>m>>itermax>>tolerance;
    INo.close();
    ifstream IN1;IN1.open("PushoverAnalysisNonlinear8SpringsDC2dof-spring.txt");
	while(! IN1.eof())
	{IN1>>DATA[X];X++;}
	IN1.close();
	N=X;
	    double F1[N/8],D1[N/8],F2[N/8],D2[N/8],F3[N/8],D3[N/8],F4[N/8],D4[N/8],f[10],Rk1[10],Rk2[10],Rk3[10],Rk4[10];
		double reaction[m],disp[m],DU[m],I[m],IT[m];
	    int Y=1,k,z,zMAX,l,Dmax;
	    double Do1,Fo1,Do2,Fo2,Do3,Fo3,Do4,Fo4;
	ifstream IN2;IN2.open("PushoverAnalysisNonlinear8SpringsDC2dof-spring.txt");
    while(IN2 >> Do1>> Fo1 >> Do2>> Fo2 >> Do3>> Fo3 >> Do4>> Fo4)
    {
    F1[Y]=Fo1;D1[Y]=Do1;
    F2[Y]=Fo2;D2[Y]=Do2;
    F3[Y]=Fo3;D3[Y]=Do3;
    F4[Y]=Fo4;D4[Y]=Do4;
    //cout<<"D["<<Y<<"]:"<<D1[Y]<<" - F["<<Y<<"]:"<<F1[Y]<<endl;
    Y++;
	}
	//cout<<"Y :"<<Y<<endl;
	    for (i=1;i<=8;i++)
    {
    	if(D4[0] < D4[i])
    	   D4[0]=D4[i];
	}
    Dmax = D4[0]; // Max displacement
    //cout<<Dmax<<endl;
	k=(Y-1);    
    cout<<green <<"     ********************************************************************\n";    
    cout<<"     *                      >> IN THE NAME OF GOD <<                    *\n";    
    cout<<"     * Pushover Analysis of Nonlinear Springs with Displacement Control *\n";
	cout<<"     *                        UNIT: [ Free Unit]                        *\n";    
    cout<<"     *------------------------------------------------------------------*\n";    
    cout<<"     *           Program is written by Salar Delavar Ghashghaei         *\n";      
    cout<<"     *                 E-mail:salar.d.ghashghaei@gmail.com              *\n";    
    cout<<"     ********************************************************************\n\n";     
    cout<<"\n\n   =========================== Input Data =============================\n"<<endl;
    cout<<"      External force [DOF(1)]:                                 "<<Force1<<endl;
    cout<<"      External force [DOF(2)]:                                 "<<Force2<<endl;
    cout<<"      Initial External Incremantal Displacement [DOF (2)]:     "<<Dini<<endl;
    cout<<"      Maximum External Incremantal Displacement [DOF (2)]:     "<<Dmax<<endl;
    cout<<"      Maximum number of iterations:                            "<<itermax<<endl;// maximum number of iterations
    cout<<"      Specified tolerance for convergence:                     "<<tolerance<<endl;//  specified tolerance for convergence
    cout<<"      Number of calculation:                                   "<<m<<endl;
        // check wrong input
  if (N%2 != 0)
	{
	cout<< red <<"\n  ======================= Error Report =========================="<<endl;
	cout<<"\n            Please check your input file! -> [ PushoverAnalysisNonlinear8SpringsDC2dof-spring.txt ]"<<endl;
	cout<<"            Imported number of input:"<<N<<endl;
	cout<<"        *** Imported number of input value is odd but it must be even ***"<<endl;
	Sleep(40000);
    }
    for(i=1;i<=N;i++)
	{
	if (DATA[i] < 0)
	{
	cout<< red <<"\n  ======================= Error Report =========================="<<endl;
	cout<<"\n            Please check your input file! -> [ PushoverAnalysisNonlinear8SpringsDC2dof-spring.txt ]"<<endl;
	cout<<"        *** Negative layer input value is not acceptable ***"<<endl;
	Sleep(40000);
	exit(1);
	}
	ifstream IN3;IN3.open("PushoverAnalysisNonlinear8SpringsDC2dof-input.txt");
	Y=0;
	while(! IN3.eof())
	{
	IN3>>CHECK[Y];
	if(CHECK[Y] < 0)
	{
	cout<< red <<"\n  ======================= Error Report =========================="<<endl;	
	cout<<"\n            Please check your input file! -> [ PushoverAnalysisNonlinear8SpringsDC2dof-input.txt ]"<<endl;
	cout<<"        *** Negative input value is not acceptable ***"<<endl;
	Sleep(40000);
	exit(1);
	}
	Y++;
	}
	IN2.close();
	}
	    cout<< yellow <<"\n   ======================= Analysis Report ========================"<<endl;
    std::clock_t start;start = clock();
    double up,duration,Kini,Ko,K[Y],Kt[Y],ff,u,du;
    double k11,k12,k22,kp,kt11,kt12,kt22,ko,F,Fi,Fii;
    u=0; // initial guess
        for (l=1;l<=8;l++)
        {
		f[l] = 0;
        }
        	for (i=1;i<=k;i++)
	        {
             Rk1[i]=(F1[i]-0)/(D1[i]-0);
             Rk2[i]=(F2[i]-F1[i])/(D2[i]-D1[i]);
             Rk3[i]=(F3[i]-F2[i])/(D3[i]-D2[i]);
             Rk4[i]=(F4[i]-F3[i])/(D4[i]-D3[i]);
             //cout<<"Rk2:"<<Rk2[i]<<endl;
            }
    cout<<"   ------------------------------------------------"<<endl;
    cout<<"     Increment    Iteration   Force  Disp.[DOF(2)] "<<endl;
    cout<<"   ------------------------------------------------"<<endl;
    for (z=1;z<=m;z++)
    {
        up = Dini*z;// Define the applied load
        for (j=1;j<=3;j++)
        {
        if (abs(f[j])>= 0 && abs(f[j])<= F1[j])
            K[j] = Rk1[j];
        else if (abs(f[j])> F1[j] && abs(f[j])<= F2[j])
			K[j] = (F1[j]+Rk2[j]*(abs(u)-D1[j]))/abs(u);
        else if (abs(f[j])> F2[j] && abs(f[j])<= F3[j])
            K[j] = (F2[j]+Rk3[j]*(abs(u)-D2[j]))/abs(u);
        else if (abs(f[j])> F3[j] && abs(f[j])<= F4[j])
            K[j] = (F3[j]+Rk4[j]*(abs(u)-D3[j]))/abs(u);
        else 
            K[j] = 0;   
        }
		for (j=4;j<=6;j++)
        {
        if (abs(f[j])>= 0 && abs(f[j])<= F1[j])
            K[j] = Rk1[j];
        else if (abs(f[j])> F1[j] && abs(f[j])<= F2[j])
            K[j] = (F1[j]+Rk2[j]*(abs(up)-abs(u)-D1[j]))/(abs(up)-abs(u));
        else if (abs(f[j])> F2[j] && abs(f[j])<= F3[j])
            K[j] = (F2[j]+Rk3[j]*(abs(up)-abs(u)-D2[j]))/(abs(up)-abs(u));
        else if (abs(f[j])> F3[j] && abs(f[j])<= F4[j])
            K[j] = (F3[j]+Rk4[j]*(abs(up)-abs(u)-D3[j]))/(abs(up)-abs(u));
        else 
            K[j] = 0;   
        }
        for (j=7;j<=8;j++)
        {
        if (abs(f[j])>= 0 && abs(f[j])<= F1[j])
            K[j] = Rk1[j];
        else if (abs(f[j])> F1[j] && abs(f[j])<= F2[j])
            K[j] = (F1[j]+Rk2[j]*(abs(up)-D1[j]))/abs(up);
        else if  (abs(f[j])> F2[j] && abs(f[j])<= F3[j])
            K[j] = (F2[j]+Rk3[j]*(abs(up)-D2[j]))/abs(up);
        else if (abs(f[j])> F3[j] && abs(f[j])<= F4[j])
            K[j] = (F3[j]+Rk4[j]*(abs(up)-D3[j]))/abs(up);
        else 
            K[j] = 0;   
        }
        k11=K[1]+K[2]+K[3]+K[4]+K[5]+K[6];
        k12=-K[2]-K[4]-K[6];
        k22=K[4]+K[5]+K[6]+K[7]+K[8];
        Fii = k12*up;
        Kini = k11;
        Fi = Force1;F=Fi-Fii;
        it = 0; // initialize iteration count
        residual = 100; // initialize residual
        while (residual > tolerance)
        {
		for (j=1;j<=3;j++)
        {
        if (abs(f[j])>= 0 && abs(f[j])<= F1[j])
            Kt[j] = Rk1[j];
        else if (abs(f[j])> F1[j] && abs(f[j])<= F2[j])
            Kt[j] = (F1[j]+Rk2[j]*(abs(u)-D1[j]))/abs(u);
        else if  (abs(f[j])> F2[j] && abs(f[j])<= F3[j])
            Kt[j] = (F2[j]+Rk3[j]*(abs(u)-D2[j]))/abs(u);
        else if (abs(f[j])> F3[j] && abs(f[j])<= F4[j])
            Kt[j] = (F3[j]+Rk4[j]*(abs(u)-D3[j]))/abs(u);
        else 
            Kt[j] = 0;   
        }
		for (j=4;j<=6;j++)
        {
        if (abs(f[j])>= 0 && abs(f[j])<= F1[j])
            Kt[j] = Rk1[j];
        else if (abs(f[j])> F1[j] && abs(f[j])<= F2[j])
            Kt[j] = (F1[j]+Rk2[j]*(abs(up)-abs(u)-D1[j]))/(abs(up)-abs(u));
        else if  (abs(f[j])> F2[j] && abs(f[j])<= F3[j])
            Kt[j] = (F2[j]+Rk3[j]*(abs(up)-abs(u)-D2[j]))/(abs(up)-abs(u));
        else if (abs(f[j])> F3[j] && abs(f[j])<= F4[j])
            Kt[j] = (F3[j]+Rk4[j]*(abs(up)-abs(u)-D3[j]))/(abs(up)-abs(u));
        else 
            Kt[j] = 0;   
        }
        for (j=7;j<=8;j++)
        {
        if (abs(f[j])>= 0 && abs(f[j])<= F1[j])
            Kt[j] = Rk1[j];
        else if (abs(f[j])> F1[j] && abs(f[j])<= F2[j])
            Kt[j] = (F1[j]+Rk2[j]*(abs(up)-D1[j]))/abs(up);
        else if  (abs(f[j])> F2[j] && abs(f[j])<= F3[j])
            Kt[j] = (F2[j]+Rk3[j]*(abs(up)-D2[j]))/abs(up);
        else if (abs(f[j])> F3[j] && abs(f[j])<= F4[j])
            Kt[j] = (F3[j]+Rk4[j]*(abs(up)-D3[j]))/abs(up);
        else 
            Kt[j] = 0;   
        }
        kt11=Kt[1]+Kt[2]+Kt[3]+Kt[4]+Kt[5]+Kt[6];
        kt12=-Kt[2]-Kt[4]-Kt[6];
        kt22=Kt[4]+Kt[5]+Kt[6]+Kt[7]+Kt[8];
        ff = kt11*u-F;
        du = (1/Kini) *(-ff); //calculate du
        u = u+du;
        residual = abs(du); // evaluate residual
        it = it + 1; // increment iteration count
        if (it == itermax)
        {
          cout<<'\t'<<z<<"\t       "<<it<<"     "<<up<<" ->   ## The solution for this step is not converged ##"<<endl;
          break;
		}
	} // while
              
        // Internal element force           
        // Stiffness Matrix for each element
        for (l=1;l<=3;l++)
        {
		f[l] = K[l]*u;
        }
        for (l=4;l<=6;l++)
        {
		f[l] = K[l]*(up-u);
        }
        for (l=7;l<=8;l++)
        {
		f[l] = K[l]*up;
        }
        // Force and Dispalcement for each increment
        reaction[z] = f[1]+f[2]+f[3]+f[7]+f[8];
		disp[z] = up;
		//disp1[z] = u;
        if (it < itermax)// iteration control
        {
        if (z < 10)
	    cout<<"\t";
        if (z >= 10 && z <= 99)
	    {cout<<'\t';cout<<'\b';}
	    if (z >= 100 && z <= 999)
	    {cout<<'\t';cout<<'\b';cout<<'\b';}	
	    cout<<z<<"\t       "<<it<<"\t"<<reaction[z]<<"     "<<up<<endl;
		}
        
        DU[z]=residual;I[z]=z;IT[z]=it;
        zMAX = z;
        if (abs(up) >= Dmax)
        {
		cout<<"\n  ## Displacement reached to ultimate displacement ##\n\n";
		break;
	    }
    } // for
    // Section bilinear fitting
    double SC,OF,AaSUM=0,UUMAX=0,FFMAX=0,k0,diy,Fy,hh[16500],Aa[zMAX+1],UU[zMAX+1],FF[zMAX+1];
    
    UU[1] = 0;FF[1] = 0;
    for (i=1;i<=zMAX;i++)
    {
	UU[i+1] = abs(disp[i]);
    FF[i+1] = abs(reaction[i]);
	}
    for (i=1;i<=zMAX;i++)
    {hh[i] = UU[i+1]-UU[i];
	 Aa[i] = (FF[i]+FF[i+1])*0.5*hh[i];
	AaSUM = AaSUM+Aa[i];
    }
	for (i=1;i<=zMAX+1;i++) // find max
	{     
	if(UUMAX < UU[i])
	{UUMAX = UU[i];}
	if(FFMAX < FF[i])
	{FFMAX = FF[i];}
	}
	k0 =FF[4] / UU[4];
    diy = (FF[i-1]*UUMAX*0.5-AaSUM)/(FF[i-1]*0.5 - k0*UUMAX*0.5);
    Fy = k0*diy;
    //cout<<diy;
    SC=UUMAX/diy;OF=FFMAX/Fy;
    cout<<"\n     Structure ductility ratio: "<<SC<<"  -   Structure over strength factor: "<<OF<<endl<<endl;
    cout<<"     ================================"<<endl;
    cout<<"     =    Bilinear curve fitted     ="<<endl;
    cout<<"     =  Displacaement    Reaction   ="<<endl;
    cout<<"     ================================"<<endl;
    cout<<"              0            0         "<<endl;
    cout<<"\t   "<<setprecision(4)<<fixed<<diy<<"\t"<<setprecision(3)<<fixed<<Fy<<endl;
    cout<<"\t   "<<setprecision(4)<<fixed<<UUMAX<<"\t"<<setprecision(3)<<fixed<<FF[i-1]<<endl;
    cout<<"     ================================\n"<<endl;
    
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	std::cout<< white <<"\n   CPU time : "<< duration <<" seconds"<<'\n';
    // current date/time based on current system
    time_t now = time(0);
    // convert now to string form
    char* dt = ctime(&now);
    cout<<"   Date and Time : " << dt;
    
    	// OUTPUT
	ofstream OUTo;OUTo.open("PushoverAnalysisNonlinear8SpringsDC2dof-output.txt");
	OUTo<<"     Date and Time : " << dt;
    OUTo<<"     ********************************************************************\n";    
    OUTo<<"     *                      >> IN THE NAME OF GOD <<                    *\n";    
    OUTo<<"     * Pushover Analysis of Nonlinear Springs with Displacement Control *\n";
	OUTo<<"     *                        UNIT: [ Free Unit]                        *\n";    
    OUTo<<"     *------------------------------------------------------------------*\n";    
    OUTo<<"     *           Program is written by Salar Delavar Ghashghaei         *\n";      
    OUTo<<"     *                 E-mail:salar.d.ghashghaei@gmail.com              *\n";    
    OUTo<<"     ********************************************************************\n\n";
    OUTo<<"   ------------------------------------------------"<<endl;
    OUTo<<"     Increment    Iteration   Force  Disp.[DOF(2)] "<<endl;
    OUTo<<"   ------------------------------------------------"<<endl;
    for (z=1;z<=zMAX;z++)
    {
    OUTo<<'\t'<<I[z]<<"\t       "<<IT[z]<<"\t"<<reaction[z]<<"     "<<disp[z]<<endl;	
    }
    cout<<"\n  - Output data is written in file -"<<endl;
    OUTo.close();
	getch();
	}
