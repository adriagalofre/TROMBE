function sol=iterative_solver_Pool(m_in,m_evap,m_out,m_ext,m_w_in,m_fr_in,cp,cpw,T,T_amb,T_cel,T_heat,T_fr,W_amb,W_evap,h_evap,beta,A_1,A_2,A_3,A_4,alpha_1_amb,alpha_1_int,alpha_2_int,alpha_3_int, alpha_4_int,epsilon,epsilon_1_star,sigma,R,F)
  
    #INICIALITZACIO DE VARIABLES
    T_w0=[293.15 293.15 293.15 293.15];                 #INITIAL GUESS
    T_w=T_w0;
    delta=1e-08;                       #TOLERANCIA SOLUCIO ITERATIVA
    err=1;                             #VAR. AUX.
    n=1;                               #VAR. AUX. 
    
	T_des=27+273.15;                   #TEMPERATURA AIRE DESHUMECTADORA. EN REALITAT DEPEN DE T_int i dels cabals
    
	#LOOP ITERACIONS SOLUCIO
    do
      #TEMPERATURA INTERIOR  
      T_int=(m_in*cp*T_des+m_evap*h_evap+alpha_1_int*T_w(1)*A_1+alpha_2_int*T_w(2)*A_2...
      +alpha_3_int*T_w(3)*A_3+alpha_4_int*T_w(4)*A_4)...
      /(m_out*cp+alpha_1_int*A_1+alpha_2_int*A_2+alpha_3_int*A_3+alpha_4_int*A_4);
      #RADIACIO TERMICA
      [qrad,g]=q_radiacio_Pool(F,epsilon,sigma,T_w);
      #RADIACIO SOLAR
      [qrads,gs]=q_radiacio_solar_Pool(F,epsilon,R,sigma,T,T_w);
      qrads_1_star=(R-1)*800+T*gs(1);
	  
      
      #TEMPERATURES PARETS
      T_w(1)=(epsilon(1)*g(1)-qrads(1)+epsilon_1_star*sigma*(T_cel^4)-qrads_1_star+alpha_1_int*T_int+alpha_1_amb*T_amb)...
              /(epsilon(1)*sigma*(T_w(1)^3)+epsilon_1_star*sigma*(T_w(1)^3)+alpha_1_int+alpha_1_amb);
      T_w(2)=(m_w_in*cpw*T_heat+m_fr_in*cpw*T_fr+epsilon(2)*g(2)-qrads(2)+alpha_2_int*T_int+m_evap*h_evap)...
              /(epsilon(2)*sigma*(T_w(2)^3)+alpha_2_int+m_w_in*cpw+m_fr_in*cpw);
			  
      T_w(3)=(epsilon(3)*g(3)-qrads(3)+alpha_3_int*T_int)...
              /(epsilon(3)*sigma*(T_w(3)^3)+alpha_3_int);
      
      T_w(4)=(epsilon(4)*g(4)-qrads(4)+alpha_4_int*T_int)...
              /(epsilon(4)*sigma*(T_w(4)^3)+alpha_4_int);
      #BALANC HUMITAT
	  #CALCUL HUMITAT
	  
	  W_des=((m_ext*W_amb)+(beta*m_out*W_int))/m_in; 
	  W_int=(m_in*W_des+m_evap*W_evap)/(m_in+m_evap);
	  
	  

	  #	CALCUL ERROR
      T_wn(n,1)=T_w(1);
      T_wn(n,2)=T_w(2);
      T_wn(n,3)=T_w(3);
      T_wn(n,4)=T_w(4);

      if n>1
      err(n)=max(abs([T_wn(n,:)-T_wn(n-1,:)]));
      else
      err(n)=max(abs((T_wn(n,:)-T_w0(:)')));
      end
      sol(n,:)=[n T_int-273.15 qrad(1) qrad(2) qrad(3) qrad(4) qrad(1)*A_1+qrad(2)*A_2+qrad(3)*A_3+qrad(4)*A_4 T_w(1)-273.15 T_w(2)-273.15 T_w(3)-273.15 T_w(4)-273.15 err(n) W_int W_des];

      if n>50
        err(n)=-1;
      end
	  
      n=n+1;
    until (err(n-1)<delta)
 end