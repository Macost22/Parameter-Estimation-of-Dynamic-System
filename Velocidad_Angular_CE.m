function [vxyz] = Velocidad_Angular_CE(Alpha,Beta,Gamma,Tiempo,Ld)
%%Derivadas en  termino de las derivadas de los angulo cardan utilizando la
%eq. 281 del libro Research Methods in biomechanic
AlphaP=zeros(3,Ld); BetaP=zeros(3,Ld); GammaP=zeros(3,Ld);
Rx=zeros(3,3,Ld);Ry=zeros(3,3,Ld);vxyz=zeros(3,Ld);
    for i=1:Ld-5
        AlphaP(:,i+1)=[(Alpha(i+2)-Alpha(i))/(Tiempo(i+2)-Tiempo(i));0;0];
        BetaP(:,i+1)=[0; (Beta(i+2)-Beta(i))/(Tiempo(i+2)-Tiempo(i)); 0];
        %GammaP(:,i+1)=[0; 0; (Gamma(i+2)-Gamma(i))/(Tiempo(i+2)-Tiempo(i))];    
        Rx(:,:,i+1)=[1 0 0; 
                    0 cos(Alpha(i+1)) sin(Alpha(i+1));
                    0 -sin(Alpha(i+1)) cos(Alpha(i+1))];

        Ry(:,:,i+1)=[cos(Beta(i+1)) 0 -sin(Beta(i+1)); 
                    0 1 0; 
                    sin(Beta(i+1)) 0 cos(Beta(i+1))]; 

        vxyz(:,i+1)=AlphaP(:,i+1)+(Rx(:,:,i+1)')*BetaP(:,i+1);%+(Rx(:,:,i+1)')*(Ry(:,:,i+1)')*GammaP(:,i+1);
    end
end

