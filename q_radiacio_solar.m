function [qrads, gs] = q_radiacio_solar(F,epsilon,R,sigma,T,T_w)
	#aux=(-1)*[0 0.7 0.7];
    #aux(1)=R;
    #aux_1= [aux(1)*F(1,:); aux(2)*F(2,:); aux(3)*F(3,:)];
    #A=eye(length(T_w))-aux_1
    A=[1-R*F(1,1) -R*F(1,2) -R*F(1,3);
		0.7*F(2,1) -1+0.7*F(2,2) 0.7*F(2,3);
		0.7*F(3,1) 0.7*F(3,2) -1+0.7*F(3,3)];
	b=[T*800 0 0];
    js=A\b';
   	gs=[F(1,:)*js F(2,:)*js F(3,:)*js];
    
	qrads=js-gs';
end