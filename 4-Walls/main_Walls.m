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

L_4=3;                                #Paret 
L_2=2;                                #Terra - Tot piscina
L_3=2;                                #Sostre - HIPOTESI ADIABATIC
L_1=sqrt((L_4^2)+(L_2-L_3)^2);        #Paret on incideix el sol (tot vidre)
W=4;                                  #Eix Z - HIPOTSESI EFECTES DESPRECIABLES
A_1=L_1*W;
A_2=L_2*W;
A_3=L_3*W;
A_4=L_4*W;
e_4=0.1;



#DADES - PROP TERMIQUES
alpha_1_int=10;
alpha_2_int=10;
alpha_3_int=10;
alpha_4_int=10;
alpha_1_amb=15;
alpha_4_hab=5;

#DADES - PROP RADIACIO
epsilon=[0.9 0.9 0.9 0.9];
epsilon_1_star=0.9;
sigma= 5.670373e-08;

#VIDRE
R=0.1;
T=0.85;
A=0.05;

#PARETS
rho_3s=0.7;
rho_4s=0.7;
lambda_4=0.5;
rho_a=1;
cp=1000;

#CONDICIONS DE CONTORN
T_hab=22+273.15;
T_amb=10+273.15;
T_cel=-10.145+273.15;
m=7*((L_4*L_3)*W+((L_2-L_3)*L_4)*W/2)*rho_a*(1/3600);


#BALANC ENERGETIC INTERIOR RECINTE
#FACTORS DE VISIO
[F]=factors_visio_Walls(L_1,L_2,L_3,L_4,A_1,A_2,A_3,A_4);


# ITERACIONS SOLUCIO TEMPERATURA
sol=iterative_solver_Walls(m,cp,T,T_amb,T_cel,T_hab,A_1,A_2,A_3,A_4,e_4,alpha_1_amb,alpha_1_int, alpha_2_int,alpha_3_int, alpha_4_int,alpha_4_hab,epsilon,epsilon_1_star,sigma,lambda_4,R,F);
				
plot(sol(:,1),sol(:,2))
title("T_{int} Convergence");
xlabel("N Iterations");
ylabel("T_{int} (C)");