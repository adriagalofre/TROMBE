function [F] = factors_visio_Walls(L_1,L_2,L_3,L_4,A_1,A_2,A_3,A_4)
	d1=sqrt((L_4^2)+(L_2^2));              #Longitudes diagonals de la seccio del volum de control
	d2=sqrt((L_4^2)+(L_3^2));
	
	#F1
	F11=0;
	F12=((L_1+L_2)-d2)/(2*L_1);
	F13=((L_1+L_3)-d1)/(2*L_1);
	F14=((d1+d2)-(L_3+L_2))/(2*L_1);
	#F2
    F21=((L_2+L_1)-d2)/(2*L_2);
	F22=0;
	F23=((d1+d2)-(L_4+L_1))/(2*L_2);
	F24=((L_2+L_4)-d1)/(2*L_2);
	#F3	
	if L_3==0
		F31=0;
		F32=0;
		F33=0;
		F34=0;
	else	
		F31=((L_3+L_1)-d1)/(2*L_3);
		F32=((d1+d2)-(L_4+L_1))/(2*L_3);
		F33=0;
		F34=((L_3+L_4)-d2)/(2*L_3);
	end
	#F4
	F41=((d1+d2)-(L_3+L_2))/(2*L_4);
	F42=((L_4+L_2)-d1)/(2*L_4);
	F43=((L_4+L_3)-d2)/(2*L_4);
	F44=0;

    F=[F11 F12 F13 F14; F21 F22 F23 F24; F31 F32 F33 F34; F41 F42 F43 F44];
    
end