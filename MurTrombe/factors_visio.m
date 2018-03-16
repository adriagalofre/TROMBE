function [F] = factors_visio(L_1,L_2,L_3,A_1,A_2,A_3)
	  F11=0;
	  F12=((L_1+L_2)-L_3)/(2*L_1);
    F13=1-F11-F12;
    F21=(A_1/A_2)*F12;
    F22=0;
    F23=1-F21-F22;
    F31=(A_1/A_3)*F13;
    F32=(A_2/A_3)*F23;
    F33=0;
    F=[F11 F12 F13; F21 F22 F23; F31 F32 F33];
    
	end