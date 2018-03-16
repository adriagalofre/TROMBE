function [qrad,g] = q_radiacio(F,epsilon,sigma,T_w)
	#aux=(1-epsilon)*(-1);
    #aux_1= [aux(1)*F(1,:); aux(2)*F(2,:); aux(3)*F(3,:)];
    #A=eye(length(T_w))-aux_1;
    A=[1-(1-epsilon(1))*F(1,1) -(1-epsilon(1))*F(1,2) -(1-epsilon(1))*F(1,3);
		-(1-epsilon(2))*F(2,1) 1-(1-epsilon(2))*F(2,2) -(1-epsilon(2))*F(2,3);
		-(1-epsilon(3))*F(3,1) -(1-epsilon(3))*F(3,2) 1-(1-epsilon(3))*F(3,3)];	
	b=sigma*[epsilon(1)*T_w(1)^4 epsilon(2)*T_w(2)^4 epsilon(3)*T_w(3)^4];
    j=A\b';
    g=[F(1,:)*j F(2,:)*j F(3,:)*j];
    #qrad=j-g';
	qrad(1)=epsilon(1)*sigma*T_w(1)^4-epsilon(1)*g(1);
	qrad(2)=epsilon(2)*sigma*T_w(2)^4-epsilon(2)*g(2);
	qrad(3)=epsilon(3)*sigma*T_w(3)^4-epsilon(3)*g(3);
	end