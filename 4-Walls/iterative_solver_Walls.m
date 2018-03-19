function sol=iterative_solver_Walls(m,cp,T,T_amb,T_cel,T_hab,A_1,A_2,A_3,A_4,e_4,alpha_1_amb,alpha_1_int, alpha_2_int,alpha_3_int, alpha_4_int,alpha_4_hab,epsilon,epsilon_1_star,sigma,lambda_4,R,F)

    #INICIALITZACIO DE VARIABLES
	if A_3==0
		T_w0=[283.15 283.15 273.15 283.15];       #INITIAL GUESS
	else
		T_w0=[283.15 283.15 283.15 283.15];       #INITIAL GUESS
	end
    
    T_w=T_w0;
    delta=1e-08;                       		  #TOLERANCIA SOLUCIO ITERATIVA
    err=1;                              	  #VAR. AUX.
    n=1;                               		  #VAR. AUX. 
    
    
	#LOOP ITERACIONS SOLUCIO
    do
      #TEMPERATURA INTERIOR  
      T_int=(m*cp*T_hab+alpha_1_int*T_w(1)*A_1+alpha_2_int*T_w(2)*A_2...
      		+alpha_3_int*T_w(3)*A_3+alpha_4_int*T_w(4)*A_4)...
      		/(m*cp+alpha_1_int*A_1+alpha_2_int*A_2+alpha_3_int*A_3+alpha_4_int*A_4);
	  
      #RADIACIO TERMICA
      [qrad,g]=q_radiacio_Walls(F,epsilon,sigma,T_w);
      #RADIACIO SOLAR
      [qrads,gs]=q_radiacio_solar_Walls(F,epsilon,R,sigma,T,T_w);
      qrads_1_star=(R-1)*800+T*gs(1);
      
      #TEMPERATURES PARETS
      T_w(1)=(epsilon(1)*g(1)-qrads(1)+epsilon_1_star*sigma*(T_cel^4)-qrads_1_star+alpha_1_int*T_int+alpha_1_amb*T_amb)...
              /(epsilon(1)*sigma*(T_w(1)^3)+epsilon_1_star*sigma*(T_w(1)^3)+alpha_1_int+alpha_1_amb);
      T_w(2)=(epsilon(2)*g(2)-qrads(2)+alpha_2_int*T_int)...
              /(epsilon(2)*sigma*(T_w(2)^3)+alpha_2_int);
      if g(3)==0
	  	T_w(3)=273.15;
	  else
	  T_w(3)=(epsilon(3)*g(3)-qrads(3)+alpha_3_int*T_int)...
              /(epsilon(3)*sigma*(T_w(3)^3)+alpha_3_int);
      end
	  T_w(4)=((T_hab/((1/alpha_4_hab)+(e_4/lambda_4)))+epsilon(4)*g(4)-qrads(4)+alpha_4_int*T_int)...
	  		/((1/((1/alpha_4_hab)+(e_4/lambda_4)))+epsilon(4)*sigma*(T_w(4)^3)+alpha_4_int);
	  
      qcond_4=(T_hab-T_w(4))/((1/alpha_4_hab)+(e_4/lambda_4));
	  
      T_wn(n,1)=T_w(1);
      T_wn(n,2)=T_w(2);
      T_wn(n,3)=T_w(3);
      T_wn(n,4)=T_w(4);
	 
      if n>1
      err(n)=max(abs([T_wn(n,:)-T_wn(n-1,:)]));
      else
      err(n)=max(abs((T_wn(n,:)-T_w0(:)')));
      end
      sol(n,:)=[n T_int-273.15 qrad(1) qrad(2) qrad(3) qrad(4) qrad(1)*A_1+qrad(2)*A_2+qrad(3)*A_3+qrad(4)*A_4 T_w(1)-273.15 T_w(2)-273.15 T_w(3)-273.15 T_w(4)-273.15 err(n)];

      if n>50
        err(n)=-1;
      end
	  
      n=n+1;
    until (err(n-1)<delta | n==50)
 end