function [axyz] = Aceleracion_Angular_CE(vxyz,Tiempo,Ld)
%Calculo de la aceleración mediante diferencias finitas centradas
axyz=zeros(3,Ld);
    for i=1:Ld-5
        axyz(:,i+1)=[(vxyz(1,i+2)-vxyz(1,i))/(Tiempo(i+2)-Tiempo(i));
                (vxyz(2,i+2)-vxyz(2,i))/(Tiempo(i+2)-Tiempo(i));
                (vxyz(3,i+2)-vxyz(3,i))/(Tiempo(i+2)-Tiempo(i))];
    end
end

