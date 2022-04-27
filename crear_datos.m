function [Datos_io, M_max,M_min] = crear_datos(genero,Alpha, Beta, axyz, L, PC)
%Cálculo de la dinámica inversa de la muñeca
Ld=length(Alpha);
if genero==true
    %distancia de la muñeca al centro de masa de la mano en cm 
    r_CD=0.3691*L-0.2; 
    %masa de la mano en kg
    rng(0,'twister');
    Ps=70;
    Pi=40;
    P=(Ps-Pi).*rand(Ld,1)+ Pi;
    m_2=zeros(Ld,1);
    for k=1:Ld
        m_2(k)=0.0061*P(k); 
    end
    m=0.0061*PC;
    %constantes para el calculo de la inercia
    b=0.285;
    c=0.182;
    a=0.233;
    %matriz de rigidez pasiva de la muñeca en Nm/kg    
    random=rand(Ld,1);
    K11i = 1.42-0.40;
    K11s = 1.42+0.40;
    K11 = (K11s-K11i).*random + K11i;
    
    K12i = -0.19-0.10;
    K12s = -0.19+0.10;
    K12 = (K12s-K12i).*random + K12i;
    
    K22i = 1.85-0.31;
    K22s = 1.85+0.31;
    K22 = (K22s-K22i).*random + K22i;
    K=zeros(2,2,Ld);
    
    for k=1:Ld
        K(:,:,k)=[K11(k) K12(k); K12(k) K22(k)];
    end
    K_min=[min(K11) min(K12); min(K12) min(K22)];
    K_max=[max(K11) max(K12); max(K12) max(K22)];
else
    r_CD=0.3502*L-0.2;
    rng(0,'twister');
    dP=20;
    Pp=67.9;
    P=dP.*randn(Ld,1) + Pp;
    m_2=zeros(Ld,1);
    for k=1:Ld
        m_2(k)=0.0056*P(k); 
    end
    m=0.0056*PC;
    b=0.241;
    c=0.152;
    a=0.206; 
    %matriz de rigidez pasiva de la muñeca en Nm/kg
    rng(0,'twister');
    random=rand(Ld,1);
    K11i = 0.96-0.33;
    K11s = 0.96+0.33;
    K11 = (K11s-K11i).*random + K11i;
    
    K12i = -0.16-0.02;
    K12s = -0.16+0.02;
    K12 = (K12s-K12i).*random + K12i;
    
    K22i = 1.49-0.32;
    K22s = 1.49+0.32;
    K22 = (K22s-K22i).*random + K22i;
    K=zeros(2,2,Ld);
    for k=1:Ld
        K(:,:,k)=[K11(k) K12(k); K12(k) K22(k)];
    end
    K_min=[min(K11) min(K12); min(K12) min(K22)];
    K_max=[max(K11) max(K12); max(K12) max(K22)];
end

M=zeros(2,Ld,Ld);M_min=zeros(2,Ld);M_max=zeros(2,Ld);
I=m*[(a*(L/100))^2 0 0; 0 (b*(L/100))^2 0; 0 0 (c*(L/100))^2]+m*[(r_CD/100)^2 0 0; 0 (r_CD/100)^2 0; 0 0 0];
g=9.776;
I_2C=zeros(3,3,Ld);
Datos_io=zeros((Ld*Ld),10);
i_aux=1;
    for i=1:Ld
        M_max(:,i)=[I(3,3) 0; 0 I(1,1)]*axyz(1:2,i)+K_max*[Alpha(i);Beta(i)]+[0; m*g*(r_CD/100)];
        M_min(:,i)=[I(3,3) 0; 0 I(1,1)]*axyz(1:2,i)+K_min*[Alpha(i);Beta(i)]+[0; m*g*(r_CD/100)];
     
        for i2=1:Ld
        I_2C(:,:,i2)=m_2(i2)*[(a*(L/100))^2 0 0; 0 (b*(L/100))^2 0; 0 0 (c*(L/100))^2]+m_2(i2)*[(r_CD/100)^2 0 0; 0 (r_CD/100)^2 0; 0 0 0];
       
        M(:,i,i2)=[I_2C(3,3,i2) 0; 0 I_2C(1,1,i2)]*axyz(1:2,i)+K(:,:,i)*[Alpha(i);Beta(i)]+[0; m_2(i2)*g*(r_CD/100)];
        Datos_io(i_aux,1)=M(1,i,i2);
        Datos_io(i_aux,2)=M(2,i,i2);
        Datos_io(i_aux,3)=Alpha(i);
        Datos_io(i_aux,4)=Beta(i);
        Datos_io(i_aux,5)=axyz(1,i);
        Datos_io(i_aux,6)=axyz(2,i);
        Datos_io(i_aux,7)=m_2(i2);
        Datos_io(i_aux,8)=K11(i);
        Datos_io(i_aux,9)=K12(i);
        Datos_io(i_aux,10)=K22(i);
        i_aux=i_aux+1;
        end
    end

end

