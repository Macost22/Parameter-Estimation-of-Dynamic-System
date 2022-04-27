function [M,M_alpha, M_beta] = Dinamica_Inversa(genero,Alpha, Beta, vxyz,axyz, L, PC)
%Cálculo de la dinámica inversa de la muñeca
Ld=length(Alpha);
if genero==true
    %distancia de la muñeca al centro de masa de la mano en cm 
    r_CD=0.3691*L-0.2; 
    %masa de la mano en kg
    m_2=0.0061*PC; 
    %constantes para el calculo de la inercia
    b=0.285;
    c=0.182;
    a=0.233;
    %matriz de rigidez pasiva de la muñeca en Nm/kg
    K=[1.42 -0.19; -0.19 1.85];%-[0.40 0.10; 0.10 0.31];   
   
else
    r_CD=0.3502*L-0.2; 
    m_2=0.0056*PC; 
    b=0.241;
    c=0.152;
    a=0.206; 
    %matriz de rigidez pasiva de la muñeca en Nm/kg
    K=[0.96 -0.16; -0.16 1.49];%-[0.33 0.02; 0.02 0.32];

end

M=zeros(2,Ld);
I_2C=m_2*[(a*(L/100))^2 0 0; 0 (b*(L/100))^2 0; 0 0 (c*(L/100))^2]+m_2*[(r_CD/100)^2 0 0; 0 (r_CD/100)^2 0; 0 0 0];

g=9.776;
M_alpha=zeros(1,Ld);M_beta=zeros(1,Ld);
r_AC=0.004;

    for i=1:Ld      
  
        M(:,i)=[I_2C(3,3) 0; 0 I_2C(1,1)]*axyz(1:2,i)+K(:,:)*[Alpha(i);Beta(i)]+[0; m_2*g*(r_CD/100)];
        M_alpha(i)=axyz(1,i)*(m_2*r_AC*(r_AC+2*(r_CD/100)*cos(Beta(i)))+I_2C(3,3)*sin(Beta(i))^2+I_2C(1,1)*cos(Beta(i))^2)+vxyz(1,i)*vxyz(2,i)*(2*(I_2C(3,3)-I_2C(1,1))*cos(Beta(i))*sin(Beta(i))-2*m_2*r_AC*(r_CD/100)*sin(Beta(i)))+Alpha(i)*K(1,1)+Beta(i)*K(1,2);
        M_beta(i)=I_2C(2,2)*axyz(2,i)-vxyz(1,i).^2*((I_2C(3,3)-I_2C(1,1))*cos(Beta(i))*sin(Beta(i))-m_2*r_AC*(r_CD/100)*sin(Beta(i)))+Alpha(i)*K(1,2)+Beta(i)*K(2,2)+m_2*g*(r_CD/100)*cos(Beta(i));
       
    end

end

