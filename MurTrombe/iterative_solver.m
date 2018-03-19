function sol=iterative_solver(m,cp,T,T_amb,T_cel, T_hab, A_1,A_2,A_3,e_2,alpha_1_amb,alpha_1_int, alpha_2_hab,alpha_2_int,alpha_3_int,epsilon,epsilon_1_star,sigma,lambda_2,R,F)
    
	#INICIALITZACIO DE VARIABLES
	T_w0=[283.15 283.15 283.15];              #INITIAL GUESS
	T_w=T_w0;
    err=1;                              	  #VAR. AUX.
    n=1;                               		  #VAR. AUX.
    delta=1e-08;                       		  #TOLERANCIA SOLUCIO ITERATIVA
    
	#LOOP ITERACIONS SOLUCIO
	do 
	  #TEMPERATURA INTERIOR	
      T_int=(m*cp*T_hab+alpha_1_int*T_w(1)*A_1+alpha_2_int*T_w(2)*A_2+alpha_3_int*T_w(3)*A_3)...
	  		/(m*cp+alpha_1_int*A_1+alpha_2_int*A_2+alpha_3_int*A_3);
 
 	  #RADIACIO TERMICA
      [qrad,g]=q_radiacio(F,epsilon,sigma,T_w);
       #RADIACIO SOLAR
	  [qrads,gs]=q_radiacio_solar(F,epsilon,R,sigma,T,T_w);
      qrads_1_star=(R-1)*800+T*gs(1);
	  
	 
      #TEMPERATURES PARETS
      T_w(1)=(epsilon(1)*g(1)-qrads(1)+epsilon_1_star*sigma*(T_cel^4)-qrads_1_star+alpha_1_int*T_int+alpha_1_amb*T_amb)...
	  		/(epsilon(1)*sigma*(T_w(1)^3)+epsilon_1_star*sigma*(T_w(1)^3)+alpha_1_int+alpha_1_amb);

      T_w(2)=((T_hab/((1/alpha_2_hab)+(e_2/lambda_2)))+epsilon(2)*g(2)-qrads(2)+alpha_2_int*T_int)...
	  		/((1/((1/alpha_2_hab)+(e_2/lambda_2)))+epsilon(2)*sigma*(T_w(2)^3)+alpha_2_int);

      T_w(3)=(epsilon(3)*g(3)-qrads(3)+alpha_3_int*T_int)...
	  		/(epsilon(3)*sigma*(T_w(3)^3)+alpha_3_int);
			
      qcond_2=(T_hab-T_w(2))/((1/alpha_2_hab)+(e_2/lambda_2));
	  
      T_wn(n,1)=T_w(1);
      T_wn(n,2)=T_w(2);
      T_wn(n,3)=T_w(3);
	  
      if n>1
      err(n)=max(abs([T_wn(n,:)-T_wn(n-1,:)]));
      else
      err(n)=max(abs((T_wn(n,:)-T_w0(:)')));
      end
      sol(n,:)=[n T_int-273.15 qrad(1) qrad(2) qrad(3) qrad(1)*A_1+qrad(2)*A_2+qrad(3)*A_3 qcond_2 T_w(1)-273.15 T_w(2)-273.15 T_w(3)-273.15 err(n)];
	  
	  if n>50
        err(n)=-1;
      end
      
	  n=n+1;
    until (err(n-1)<delta | n==50)
 end