clc
clear all

#### MUR TROMBE-POOL ####
#########################

##HIPOTESIS:
#
#
#
#
#
#
#

#DADES - GEOMETRIA

L_1=10;                                #Paret on incideix el sol (tot vidre)
L_2=30;                                #Terra - Tot piscina
L_3=30;                                #Sostre - HIPOTESI ADIABATIC
L_4=10;                                #Paret  - HIPOTESI ADIABATICA
W=30;                                  #Eix Z - HIPOTSESI EFECTES DESPRECIABLES
A_1=L_1*W;
A_2=L_2*W;
A_3=L_3*W;
A_4=L_4*W;
d1=sqrt((L_1^2)+(L_3^2));              #Longitudes diagonals de la seccio del volum de control
d2=d1;


#DADES - PROP TERMIQUES
alpha_1_int=10;
alpha_2_int=5;
alpha_3_int=10;
alpha_4_int=10;
alpha_1_amb=15;

#DADES - PROP RADIACIO
epsilon=[0.9 0.2 0.9 0.9];
epsilon_1_star=0.9;
sigma= 5.670373e-08;

#AIGUA PISCINA
h_evap=2500;                           #CALOR LATENT DE VAPORITZACIO. APROXIMAT. KJ/KG

#VIDRE
R=0.1;
T=0.85;
A=0.05;

#PARETS
rho_a=1;
rho_3s=0.7;
rho_4s=0.7;
lambda_3=0.5;
rho_air=1;
cp=1000;

#CONDICIONS DE CONTORN
T_hab=22;
T_amb=10;
T_cel=-10.145;
m_in=7*(((L_2*L_3)/2)*W)*rho_a*(1/3600);
m_evap=10;                                       #HIPOTESI CABAL EVAPORACIO CONSTANT
m_out=m_in+m_evap;

#BALANÇ ENERGETIC INTERIOR RECINTE
#FACTORS DE VISIO
[F]=factors_visio_Pool(L_1,L_2,L_3,L_4,A_1,A_2,A_3,A_4,d1,d2);


# ITERACIONS SOLUCIO TEMPERATURA
sol=iterative_solver_Pool(m_in,m_evap,m_out,cp,T,T_amb,T_cel,h_evap,A_1,A_2,A_3,A_4,alpha_1_amb,alpha_1_int, alpha_2_int,alpha_3_int, alpha_4_int,epsilon,epsilon_1_star,sigma,lambda_3,R,F)
