function [qrad,g] = q_radiacio_Pool(F,epsilon,sigma,T_w)
	  aux=(1-epsilon)*(-1);
    aux_1= [aux(1)*F(1,:); aux(2)*F(2,:); aux(3)*F(3,:); aux(4)*F(4,:)];
    A=eye(length(T_w))-aux_1;
    b=sigma*[epsilon(1)*T_w(1)^4; epsilon(2)*T_w(2)^4; epsilon(3)*T_w(3)^4; epsilon(4)*T_w(4)^4];
    j=A\b;
    g=[F(1,:)*j; F(2,:)*j; F(3,:)*j; F(4,:)*j];
    qrad=j-g;
	end