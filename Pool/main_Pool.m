clc
clear all

########## MUR TROMBE-POOL ##########
#####################################

######HIPOTESIS:######
#
#
#
#
#
#
#

######DADES - GEOMETRIA######
#
#	  ___L_3__
#    |		  |
#    |        |
#L_4 |        | L_1/Vidre
#    |________| 
#     L_2/Pool
#

L_4=30;                                #Paret on incideix el sol (tot vidre)
L_2=20;                                #Terra - Tot piscina
L_3=20;                                #Sostre - HIPOTESI ADIABATIC
L_1=sqrt((L_4^2)+(L_2-L_3)^2);         #Paret  - HIPOTESI ADIABATICA
W=4;                                   #Eix Z - HIPOTSESI EFECTES DESPRECIABLES
A_1=L_1*W;
A_2=L_2*W;
A_3=L_3*W;
A_4=L_4*W;
A_pool=A_2;

######DADES - PROP TERMIQUES######
alpha_1_int=10;
alpha_2_int=10;
alpha_3_int=10;
alpha_4_int=10;
alpha_1_amb=15;

######DADES - PROP RADIACIO######
epsilon=[0.9 0.95 0.9 0.9];
epsilon_1_star=0.9;
sigma= 5.670373e-08;

######AIGUA PISCINA#######
h_evap=2500;                            #CALOR LATENT DE VAPORITZACIO. APROXIMAT. kJ/KG
#h_evap=0;
W_evap=0.9
#Bernier Formula
n=(5)/A_2;
N=0;
#m_evap=1*m_in;                         #HIPOTESI CABAL EVAPORACIO CONSTANT/
#m_evap=0;
Ga=0.65;
We=0.020;
Was=0.0225;
m_evap=(A_pool*((16+133*n)*(We-Ga*Was))+0.1*N)/3600            #Bernier Formula with: n=nadadors/m², N=espectadors, Ga=Grau saturació, We(kgag/Kga) H.Abs del aire saturat
#															   # a la temperatura de l'aigua de la piscina. Was H.Abs del aire saturat a la temperatura de l'aire interior.
#Auxiliary Heating
T_fr=10+273.15;					    #25 aigua que és busca al vas de la piscina
cpw=4.1814;        					    # Kj/Kg·K
m_fr_in=500*m_evap
#
#Auxiliary Heating
T_heat=25+273.15;					    #25 aigua que és busca al vas de la piscina
cpw=4.1814;        					    # Kj/Kg·K
m_h_in=100*m_evap
	
######VIDRE######
R=0.1;
T=0.85;
A=0.05;
#PARET
rho_3s=0.7;
rho_4s=0.7;
rho_a=1;
cp=1000;

######DESHUMECTADORA######
beta=rand();                          #Factor by pass. La funció rand() retorna un número entre 0 i 1.
m_out=7*((L_4*L_3)*W)*rho_a*(1/3600)
m_in=m_out-m_evap
m_ext=m_in-beta*m_out

######CONDICIONS DE CONTORN######
W_amb=0.027125;                       #Kgw/kga  
T_amb=30+273.15;
T_cel=-10.145+273.15;



######BALANC ENERGETIC INTERIOR RECINTE######
#
######FACTORS DE VISIO######
[F]=factors_visio_Pool(L_1,L_2,L_3,L_4,A_1,A_2,A_3,A_4);


######ITERACIONS SOLUCIO TEMPERATURA######
 
	sol=iterative_solver_Pool(m_in,m_evap,m_out,m_ext,m_h_in,m_fr_in,cp,cpw,T,T_amb,T_cel,T_heat,T_fr,W_amb,W_evap,h_evap,beta,A_1,A_2,A_3,A_4,alpha_1_amb,alpha_1_int, alpha_2_int,alpha_3_int, alpha_4_int,epsilon,epsilon_1_star,sigma,R,F);
	
	subplot(2,1,1)
	plot(sol(:,1),sol(:,2),";T_{int};",sol(:,1),sol(:,8),";T_{w1};",sol(:,1),sol(:,9),";T_{w2};",sol(:,1),sol(:,10),";T_{w3};",sol(:,1), sol(:,11),";T_{w4};");
	title("T_{int} Convergence");
	xlabel("N Iterations");
	ylabel("T (C)");
	subplot(2,1,2)
	plot(sol(:,1),sol(:,13),";W_{int};",sol(:,1), sol(:,14),";W_{des};");
	title("W_{int} Convergence");
	xlabel("N Iterations");
	ylabel("Humidity ratio (W)");
	
	#legend("T_{int}","T_{w1}","T_{w2}","T_{w3}","T_{w4}");
	#h=legend("show");