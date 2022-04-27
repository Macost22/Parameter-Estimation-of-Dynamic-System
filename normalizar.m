function [Datos_n] = normalizar(Datos,min12,max12,min34,max34,min56,max56,min7,max7,min810,max810)
%función para normalizar los datos
Datos_n=zeros(size(Datos));
Ld=length(Datos(:,1));
% min12=min(Datos(:,1:2),[],'all');
% max12=max(Datos(:,1:2),[],'all');
% 
% min34=min(Datos(:,3:4),[],'all');
% max34=max(Datos(:,3:4),[],'all');
% 
% min56=min(Datos(:,5:6),[],'all');
% max56=max(Datos(:,5:6),[],'all');
% 
% min7=min(Datos(:,7));
% max7=max(Datos(:,7));
% 
% min810=min(Datos(:,8:10),[],'all');
% max810=max(Datos(:,8:10),[],'all');
    
 for i=1:Ld   
  Datos_n(i,1)=(Datos(i,1) - min12) / ( max12 - min12); 
  Datos_n(i,2)=(Datos(i,2) - min12) / ( max12 - min12); 
  Datos_n(i,3)=(Datos(i,3) - min34) / ( max34 - min34); 
  Datos_n(i,4)=(Datos(i,4) - min34) / ( max34 - min34); 
  Datos_n(i,5)=(Datos(i,5) - min56) / ( max56 - min56); 
  Datos_n(i,6)=(Datos(i,6) - min56) / ( max56 - min56); 
  Datos_n(i,7)=(Datos(i,7) - min7) / ( max7 - min7); 
  Datos_n(i,8)=(Datos(i,8) - min810) / ( max810 - min810); 
  Datos_n(i,9)=(Datos(i,9) - min810) / ( max810 - min810); 
  Datos_n(i,10)=(Datos(i,10) - min810) / ( max810 - min810); 
 end        
   
end

