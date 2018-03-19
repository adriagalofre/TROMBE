function [qrads, gs] = q_radiacio_solar_Walls(F,epsilon,R,sigma,T,T_w)
	#aux=(1-epsilon)*(-1);
    #aux(1)=R;
    #aux_1= [aux(1)*F(1,:); aux(2)*F(2,:); aux(3)*F(3,:); aux(4)*F(4,:)];
    #A=eye(length(T_w))-aux_1;
    #A=eye(length(T_w))-aux_1
	#Assumim rho superficie-sol 0.7 
    A=[1-R*F(1,1) -R*F(1,2) -R*F(1,3) -R*F(1,4);
		0.7*F(2,1) -1+0.7*F(2,2) 0.7*F(2,3) 0.7*F(2,4);
		0.7*F(3,1) 0.7*F(3,2) -1+0.7*F(3,3) 0.7*F(3,4);
		0.7*F(4,1) 0.7*F(4,2) 0.7*F(4,3) -1+0.7*F(4,4)];
	
	b=[T*800 0 0 0];
    j=A\b';
    gs=[F(1,:)*j F(2,:)*j F(3,:)*j F(4,:)*j];
    qrads=j-gs';
end