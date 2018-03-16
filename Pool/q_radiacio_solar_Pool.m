function [qrads, gs] = q_radiacio_solar_Pool(F,epsilon,R,sigma,T,T_w)
	  aux=(1-epsilon)*(-1);
    aux(1)=R;
    aux_1= [aux(1)*F(1,:); aux(2)*F(2,:); aux(3)*F(3,:); aux(4)*F(4,:)];
    A=eye(length(T_w))-aux_1;
    b=[T*800 0 0 0];
    j=A\b';
    gs=[F(1,:)*j F(2,:)*j F(3,:)*j F(4,:)*j];
    qrads=j-gs';
	end