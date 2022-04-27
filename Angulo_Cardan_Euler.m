function [Alpha,Beta,Gamma] = Angulo_Cardan_Euler(R,Ld)
%CALCULO DE LOS ÁNGULOS DE CARDA EULER A PARTIR DE LA MATRIZ DE ROTACIÓN
Alpha=zeros(Ld,1); Beta=zeros(Ld,1); Gamma=zeros(Ld,1);
    for i=1:Ld
        Alpha(i)=atan2(-R(3,2,i),R(3,3,i));
        Beta(i)=atan2(R(3,1,i),sqrt(R(1,1,i)^2+R(2,1,i)^2));
        Gamma(i)=atan2(-R(2,1,i),R(1,1,i));
    end
end

