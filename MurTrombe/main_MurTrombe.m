clc
clear all

#### MUR TROMBE ####
####################

#DADES - GEOMETRIA
L_2=3;
L_3=2;
L_1=sqrt(L_2^2+L_3^2);
W=4;
A_1=L_1*W;
A_2=L_2*W;
A_3=L_3*W;
e_2=0.1;

#DADES - PROP TERMIQUES
alpha_1_int=10;
alpha_2_int=10;
alpha_3_int=10;
alpha_1_amb=15;
alpha_2_hab=5;

epsilon=[0.9 0.9 0.9];
epsilon_1_star=0.9;
sigma= 5.670373e-08;

#VIDRE
R=0.1;
T=0.85;
A=0.05;
#PARETS
rho_2s=0.7;
rho_3s=0.7;
lambda_2=0.5;
rho_a=1;
cp=1000;

#CONDICIONS DE CONTORN
T_hab=22+273.15;
T_amb=10+273.15;
T_cel=-10.145+273.15;
m=7*(((L_2*L_3)/2)*W)*rho_a*(1/3600);

#BALANÇ ENERGETIC INTERIOR RECINTE
#FACTORS DE VISIO
[F]=factors_visio(L_1,L_2,L_3,A_1,A_2,A_3);


# ITERACIONS SOLUCIO TEMPERATURA
sol=iterative_solver(m,cp,T,T_amb,T_cel, T_hab, A_1,A_2,A_3,e_2,alpha_1_amb,alpha_1_int, alpha_2_hab,alpha_2_int,alpha_3_int,epsilon,epsilon_1_star,sigma,lambda_2,R,F)

