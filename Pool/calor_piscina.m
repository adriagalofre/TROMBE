function [qheat] = calor_piscina(T_w(2),m_evap)
	#funcio balanc piscina + heater
	T_heat=25+273.15;					#25 aigua que és busca al vas de la piscina
	cpw=4.1814;        					# Kj/Kg·K
	m_w_in=m_evap;
	
	qheat=m_w_in*cpw*(T_w(2)-T_heat);
end